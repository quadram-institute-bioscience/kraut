from typing import List, Optional, Dict, Union
import sys

class KrakenNode:
    def __init__(
        self,
        percent: float,
        clade_counts: int,
        taxon_counts: int,
        rank_code: str,
        tax_id: int,
        name: str,
        depth: int = 0,
        parent: Optional['KrakenNode'] = None
    ):
        self.percent = percent
        self.clade_counts = clade_counts
        self.taxon_counts = taxon_counts
        self.rank_code = rank_code
        self.tax_id = tax_id
        self.name = name
        self.depth = depth
        self.parent = parent
        self.children: List['KrakenNode'] = []
    
    def add_child(self, child: 'KrakenNode'):
        child.parent = self
        self.children.append(child)

    def get_children(self) -> List['KrakenNode']:
        return self.children

    def to_string(self) -> str:
        """
        Returns the string representation of this node line, 
        attempting to reproduce Kraken2 output format.
        """
        # Format: percent (6.2f), clade_counts, taxon_counts, rank_code, tax_id, indented_name
        # Kraken uses strict tabulation.
        # Example:
        # 100.00	3	0	R	1	root
        #  74.59	2173760	2173760	U	0	unclassified
        
        # Name indentation is 2 spaces per depth
        indent = "  " * self.depth
        line = f"{self.percent:6.2f}\t{self.clade_counts}\t{self.taxon_counts}\t{self.rank_code}\t{self.tax_id}\t{indent}{self.name}"
        return line


class KrakenReport:
    def __init__(self):
        self.root: Optional[KrakenNode] = None
        self.unclassified: Optional[KrakenNode] = None
        self.nodes: Dict[int, KrakenNode] = {}
        # We might need to store order if we want to reproduce exactly, 
        # but tree traversal usually reconstructs it.
        # However, unclassified is usually first in file.

    @classmethod
    def from_file(cls, path: str) -> 'KrakenReport':
        report = cls()
        
        # We need to maintain a stack to reconstruct the tree structure based on indentation/depth
        # But Kraken report indentation is strictly rank based or depth based?
        # Actually it's just visual. The tree structure is implicit.
        # Logic:
        # "root" (taxid 1) is the root of the taxonomy.
        # "unclassified" (taxid 0) is separate.
        # Children follow parents. Depth increases.
        
        stack: List[KrakenNode] = []
        
        with open(path, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) != 6:
                    continue # Should probably error or warn
                
                percent = float(parts[0].strip())
                clade_counts = int(parts[1])
                taxon_counts = int(parts[2])
                rank_code = parts[3]
                tax_id = int(parts[4])
                raw_name = parts[5]
                
                # Determine depth from leading spaces of name
                # Kraken uses 2 spaces per level of indentation
                stripped_name = raw_name.lstrip(' ')
                num_leading_spaces = len(raw_name) - len(stripped_name)
                depth = num_leading_spaces // 2
                
                node = KrakenNode(
                    percent=percent,
                    clade_counts=clade_counts,
                    taxon_counts=taxon_counts,
                    rank_code=rank_code,
                    tax_id=tax_id,
                    name=stripped_name,
                    depth=depth
                )
                
                report.nodes[tax_id] = node
                
                if tax_id == 0 or node.name == 'unclassified':
                    report.unclassified = node
                    # Unclassified is not part of the main tree structure usually, 
                    # or it's a sibling of root?
                    # The prompt says: "Treat unclassified as a special class extra tree"
                    continue
                    
                if tax_id == 1 or node.name == 'root':
                    report.root = node
                    stack = [node]
                    continue
                
                # Tree building logic
                # Find parent in stack
                # If current depth > stack top depth, then stack top is parent
                # If current depth == stack top depth, then stack[-2] is parent (sibling)
                # If current depth < stack top depth, pop until find parent (depth - 1)
                
                if not stack:
                     # Should not happen if root comes first (after unclassified)
                     # But some files might be fragments. 
                     # For now assume valid kraken report starting with root or unclassified
                     pass

                # Adjust stack to find parent
                # We want the node at depth - 1 to be the parent
                while stack and stack[-1].depth >= depth:
                    stack.pop()
                
                if stack:
                    parent = stack[-1]
                    parent.add_child(node)
                
                stack.append(node)

        return report

    def to_string(self, 
                  min_fract: float = 0.0, 
                  min_count: int = 0, 
                  min_level: Optional[str] = None, 
                  max_level: Optional[str] = None) -> str:
        
        output_lines = []
        
        # 1. Print Unclassified if present
        # Assuming we don't filter unclassified unless specifically requested?
        # The prompt for single-report says "By default the output will be the same as the input"
        if self.unclassified:
            # Check filters for unclassified too? Usually YES.
            if self._satisfies_filter(self.unclassified, min_fract, min_count, min_level, max_level):
                 output_lines.append(self.unclassified.to_string())
        
        # 2. Print Tree
        if self.root:
             self._traverse_print(self.root, output_lines, min_fract, min_count, min_level, max_level)
             
        return "\n".join(output_lines) + "\n"

    def _traverse_print(self, 
                        node: KrakenNode, 
                        lines: List[str],
                        min_fract: float, 
                        min_count: int, 
                        min_level: Optional[str], 
                        max_level: Optional[str]):
        
        # Check filters
        # For max-level (L), we assume "do not print below this level but aggregate".
        # But this method is for printing an existing tree. 
        # Aggregation logic usually requires re-calculating counts if we prune.
        # BUT, standard Kraken reports are cumulative. 
        # If we just stop printing at max-level, the counts at that level already include children.
        # So we just prune the traversal.
        
        if not self._satisfies_filter(node, min_fract, min_count, min_level, max_level):
            return

        lines.append(node.to_string())
        
        # If we hit max_level, we stop recursing?
        # TODO: Implement robust level comparison. 
        # For now, minimal implementation.
        
        for child in node.children:
            self._traverse_print(child, lines, min_fract, min_count, min_level, max_level)

    def _satisfies_filter(self, node: KrakenNode, min_fract, min_count, min_level, max_level) -> bool:
        if node.percent < min_fract:
            return False
        if node.clade_counts < min_count:
            return False
        
        # Level filtering requires rank parsing (U, R, D, K, P, C, O, F, G, S)
        # and maybe standardizing (S1, S2 etc).
        # We'll leave strict level filtering for refinement step to ensure we get 
        # the order right.
        
        return True
