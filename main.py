import numpy as np

from sklearn.linear_model import LinearRegression


import json


tokens = {}

with open("tokens.json") as f:
    tokens = json.load(f) 


def get_token(sym: str):
    global tokens
    
    return tokens[sym]


def atom2vec(symbol: str) -> np.ndarray:
    a = tokens[symbol]

    return np.array([
        a['en'] / 4,
        a['valence'] / 8,
        a['cov_radius'] / 2,
        a['lone_pairs'] / 6
    ], dtype=np.float32)

def bond2vec(sym_a: str, sym_b: str, bond_order: int) -> np.ndarray:
    """將鍵結轉換為特徵向量（3 維）。"""


    a = tokens[sym_a]
    b = tokens[sym_b]
    
    return np.array([
        bond_order / 3,

        abs(
            a['en']
            - b['en']
        ) / 4,

        (a['cov_radius'] + b['cov_radius']) / 4,
    ], dtype=np.float32)



def mol2vec(mol: dict) -> np.ndarray:
    """將完整分子轉換成固定長度的特徵向量。"""

    atoms = mol["atoms"]
    edges = mol["edges"]

    if not atoms:
        raise ValueError("分子不能沒有原子")

    # 1. 檢查原子符號
    for symbol in atoms:
        if symbol not in tokens:
            raise ValueError(f"tokens.json 沒有元素：{symbol}")

    # 2. 原子特徵：每個原子 4 維
    atom_vecs = np.stack([
        atom2vec(symbol)
        for symbol in atoms
    ])

    atom_mean = atom_vecs.mean(axis=0)
    atom_std = atom_vecs.std(axis=0)

    # 3. 鍵結特徵：每條鍵結 3 維
    bond_vecs = []

    for a, b, order in edges:
        if not (0 <= a < len(atoms)):
            raise ValueError(f"原子索引錯誤：{a}")

        if not (0 <= b < len(atoms)):
            raise ValueError(f"原子索引錯誤：{b}")

        if a == b:
            raise ValueError("不允許原子與自己形成鍵結")

        if order not in (1, 2, 3):
            raise ValueError(f"不支援的鍵級：{order}")

        bond_vecs.append(
            bond2vec(atoms[a], atoms[b], order)
        )

    if bond_vecs:
        bond_vecs = np.stack(bond_vecs)

        bond_mean = bond_vecs.mean(axis=0)
        bond_std = bond_vecs.std(axis=0)
    else:
        bond_mean = np.zeros(3, dtype=np.float32)
        bond_std = np.zeros(3, dtype=np.float32)

    # 4. 元素組成：固定順序，每種元素一個特徵
    element_symbols = sorted(tokens.keys())

    element_counts = np.array([
        atoms.count(symbol)
        for symbol in element_symbols
    ], dtype=np.float32)

    element_fractions = element_counts / len(atoms)

    # 5. 鍵級分布：單鍵、雙鍵、三鍵的比例
    bond_order_fractions = np.zeros(3, dtype=np.float32)

    for _, _, order in edges:
        bond_order_fractions[order - 1] += 1

    if edges:
        bond_order_fractions /= len(edges)

    # 6. 分子大小
    size_features = np.array([
        np.log1p(len(atoms)),
        np.log1p(len(edges)),
    ], dtype=np.float32)

    # 7. 組合成固定長度向量
    return np.concatenate([
        atom_mean,               # 4 維
        atom_std,                # 4 維
        bond_mean,               # 3 維
        bond_std,                # 3 維
        element_fractions,       # len(tokens) 維
        bond_order_fractions,    # 3 維
        size_features,           # 2 維
    ]).astype(np.float32)

def get_neighbors(mol):
    atoms = mol["atoms"]
    neighbors = [[] for _ in atoms]

    for a, b, order in mol["edges"]:
        neighbors[a].append((atoms[b], order))
        neighbors[b].append((atoms[a], order))

    return neighbors

methane = {
    "atoms": ["C","H","H","H", "H"],
    "edges": [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
}

methanol = {
    "atoms": ["C", "H", "H", "H", "O", "H"],
    "edges": [
        (0, 1, 1),
        (0, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (4, 5, 1),
    ],
}
ethanol = {
    "atoms": ["C", "H", "H",  "H", "C", "H", "H", "O", "H"],
    "edges": [
        (0, 1, 1),
        (0, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (4, 5, 1),
        (4, 6, 1),
        (4, 7, 1),
        (7, 8, 1),
    ],
}

dimethyl_ether = {
    "atoms": ["C", "H", "H",  "H", "O", "C", "H", "H",  "H"],
    "edges": [
        (0, 1, 1),
        (0, 2, 1),
        (0, 3, 1),
        (0, 4, 1),
        (4, 5, 1),
        (5, 6, 1),
        (5, 7, 1),
        (5, 8, 1),
    ],
}

ammonia = {
    "atoms": ["N","H","H","H"],
    "edges": [(0,1,1),(0,2,1),(0,3,1)],
}
water = {
    "atoms": ["O","H","H"],
    "edges": [(0,1,1),(0,2,1)],
}

molecules = [methane, methanol, ammonia, ethanol, dimethyl_ether, water]

X = np.stack([
    mol2vec(mol)
    for mol in molecules
])

print("X 的形狀：", X.shape)

def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """計算兩向量的餘弦相似度。"""
    norm_a, norm_b = np.linalg.norm(a), np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        raise ValueError("零向量無法計算餘弦相似度")
    return float(np.dot(a, b) / (norm_a * norm_b))


print("餘弦相似度：", cosine_similarity(X[3], X[4]))

for name, mol in [
    ("ethanol", ethanol),
    ("dimethyl_ether", dimethyl_ether),
]:
    neighbors = get_neighbors(mol)

    print(f"\n{name}")

    for i, atom in enumerate(mol["atoms"]):
        if atom == "O":
            print(
                f"O at index {i}: "
                f"{neighbors[i]}"
            )

assert X.ndim == 2
assert X.shape[0] == len(molecules)
assert np.isfinite(X).all()
assert all(
    len(mol2vec(mol)) == X.shape[1]
    for mol in molecules
)