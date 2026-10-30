# 分子段落
from dataclasses import dataclass


@dataclass
class MolFragment:
    atoms: list[str]
    
    edges: list[tuple]
    
    @staticmethod
    def merge2(f1: 'MolFragment', f2: 'MolFragment') -> 'MolFragment':
        # 1. 複製 f1 的資料
        new_atoms = f1.atoms.copy()
        new_edges = f1.edges.copy()
        
        # 2. 計算 f2 的索引偏移量，並將 f2 的邊加上偏移量後併入
        offset = len(f1.atoms)
        new_atoms.extend(f2.atoms)
        
        for idx1, idx2, b_order in f2.edges:
            new_edges.append((idx1 + offset, idx2 + offset, b_order))
            
        # 3. 尋找對接點 (標記為數字的原子)
        idx_f1 = -1
        idx_f2 = -1
        
        for i, atom in enumerate(f1.atoms):
            if atom.isdigit():
                idx_f1 = i
                break
                
        for i, atom in enumerate(f2.atoms):
            if atom.isdigit():
                idx_f2 = i + offset  # 記得加上偏移量
                break
        
        # 4. 如果找到對接點，將兩邊的鄰居連接起來，並移除對接點
        if idx_f1 != -1 and idx_f2 != -1:
            # 找到與 f1 標記點連接的原子與鍵級
            nb_f1, bo_f1 = [(e[1] if e[0] == idx_f1 else e[0], e[2]) for e in new_edges if idx_f1 in (e[0], e[1])][0]
            # 找到與 f2 標記點連接的原子與鍵級
            nb_f2, bo_f2 = [(e[1] if e[0] == idx_f2 else e[0], e[2]) for e in new_edges if idx_f2 in (e[0], e[1])][0]
            
            # 移除舊的、連到標記點的邊
            new_edges = [e for e in new_edges if idx_f1 not in (e[0], e[1]) and idx_f2 not in (e[0], e[1])]
            
            # 建立新的邊連接兩個片段 (這裡假設鍵級不變，取 bo_f1)
            new_edges.append((nb_f1, nb_f2, bo_f1))
            
            # 移除原子列表中的數字標記（由後往前刪除以防索引失效）
            # 因為 idx_f2 一定大於 idx_f1
            new_atoms.pop(idx_f2)
            new_atoms.pop(idx_f1)
            
            # 由於移除了兩個原子，所有大於被刪除索引的邊索引都需要往前修正
            # 為了簡化此處的索引重排，我們可以直接重建邊的索引
            # 但更簡單的做法是：在分子圖層面，通常會把標記點視為「虛擬原子」或直接用特定的連接機制處理。
            
            # 修正因 pop 導致的邊索引錯位
            remap = {}
            curr = 0
            for i in range(len(f1.atoms) + len(f2.atoms)):
                if i == idx_f1 or i == idx_f2:
                    continue
                remap[i] = curr
                curr += 1
                
            adjusted_edges = [(remap[e[0]], remap[e[1]], e[2]) for e in new_edges]
            new_edges = adjusted_edges

        return MolFragment(new_atoms, new_edges)
    

if __name__ == "__main__":
    # 羥基
    HYDROXYL = MolFragment(["1", "O", "H"], [(0, 1, 1),(1, 2, 1)])

    # 甲基
    METHYL = MolFragment( ["C","H","H","H", "1"],  [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])

    print(MolFragment.merge2(METHYL, HYDROXYL))

