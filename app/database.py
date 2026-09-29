"""Thread-safe read-only, parameterized queries over the curated local store."""
from contextlib import closing
from pathlib import Path
import sqlite3

from .chemistry import _NAMES as _CURATED_MOLECULE_NAMES


DEFAULT_DATABASE_PATH = Path(__file__).resolve().parents[1] / "data" / "chemistry.sqlite"
TARGETS = frozenset({"logS", "logD74", "hydration_free_energy"})
# This deliberately small translation table resolves only names already in the
# chemistry module's curated offline map. It is not a general name-to-structure
# resolver, and aliases are never expanded into SMILES substring searches.
_CHINESE_NAME_ALIASES = {
    "阿司匹林": "Aspirin",
    "咖啡因": "Caffeine",
    "乙醇": "Ethanol",
    "对乙酰氨基酚": "Paracetamol",
    "苯": "Benzene",
}
_CANONICAL_BY_KNOWN_NAME = {name: smiles for smiles, name in _CURATED_MOLECULE_NAMES.items()}


class MoleculeDatabase:
    def __init__(self, path: Path | str | None = None):
        self.path = Path(path).resolve() if path is not None else DEFAULT_DATABASE_PATH.resolve()
        if not self.path.is_file():
            raise FileNotFoundError("本地分子数据库不存在；请运行 python scripts/prepare_data.py。")
        with closing(self._connect()) as connection:
            tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if not {"molecules", "observations", "datasets"}.issubset(tables):
                raise ValueError("本地分子数据库结构不匹配；请重新运行数据准备脚本。")

    def _connect(self):
        # Each request owns its connection; no shared worker-thread connection.
        connection = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True, timeout=10)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only = ON")
        return connection

    def stats(self):
        with closing(self._connect()) as connection:
            targets = {row["target"]: {"observations": row["observations"], "molecules": row["molecules"], "unit": row["unit"]}
                for row in connection.execute("SELECT target,COUNT(*) AS observations,COUNT(DISTINCT molecule_id) AS molecules,MIN(unit) AS unit FROM observations GROUP BY target ORDER BY target")}
            datasets = [dict(row) for row in connection.execute("SELECT * FROM datasets ORDER BY id")]
            return dict(molecules=connection.execute("SELECT COUNT(*) FROM molecules").fetchone()[0],
                observations=connection.execute("SELECT COUNT(*) FROM observations").fetchone()[0],
                dataset_count=len(datasets), datasets=datasets, targets=targets,
                target_counts={target: row["observations"] for target, row in targets.items()})

    def search(self, query="", target=None, limit=20, offset=0):
        if not isinstance(query, str) or len(query) > 2000:
            raise ValueError("query 必须为最多 2000 个字符的字符串。")
        if target is not None and target not in TARGETS:
            raise ValueError("未知的性质筛选项。")
        if type(limit) is not int or not 1 <= limit <= 100:
            raise ValueError("limit 必须为 1 到 100 的整数。")
        if type(offset) is not int or not 0 <= offset <= 1_000_000:
            raise ValueError("offset 必须为 0 到 1000000 的整数。")
        filters, params = [], []
        ordering, ordering_params = "m.id", []
        query = query.strip()
        if query:
            escaped = query.replace("!", "!!").replace("%", "!%").replace("_", "!_")
            pattern = f"%{escaped}%"
            search_filter = "(m.name LIKE ? ESCAPE '!' OR m.canonical_smiles LIKE ? ESCAPE '!' OR m.formula LIKE ? ESCAPE '!' OR EXISTS(SELECT 1 FROM observations AS q WHERE q.molecule_id=m.id AND q.source_id LIKE ? ESCAPE '!')"
            params.extend([pattern] * 4)
            known_canonical = _CANONICAL_BY_KNOWN_NAME.get(_CHINESE_NAME_ALIASES.get(query))
            if known_canonical is not None:
                search_filter += " OR m.canonical_smiles = ?"
                params.append(known_canonical)
            filters.append(search_filter + ")")
            ordering = "CASE WHEN m.canonical_smiles = ? THEN 0 WHEN m.name COLLATE NOCASE = ? THEN 1 ELSE 2 END, m.id"
            ordering_params = [known_canonical or query, query]
        if target is not None:
            filters.append("EXISTS(SELECT 1 FROM observations AS t WHERE t.molecule_id=m.id AND t.target=?)")
            params.append(target)
        where = " WHERE " + " AND ".join(filters) if filters else ""
        with closing(self._connect()) as connection:
            total = connection.execute("SELECT COUNT(*) FROM molecules AS m" + where, params).fetchone()[0]
            rows = connection.execute("SELECT m.id,m.name,m.canonical_smiles,m.formula,m.mw FROM molecules AS m" + where + " ORDER BY " + ordering + " LIMIT ? OFFSET ?", params + ordering_params + [limit, offset]).fetchall()
            items = []
            for row in rows:
                item = dict(row)
                item["targets"] = [value[0] for value in connection.execute("SELECT DISTINCT target FROM observations WHERE molecule_id=? ORDER BY target", (item["id"],))]
                items.append(item)
            return dict(total=total, items=items, limit=limit, offset=offset)

    def lookup(self, canonical_smiles):
        if not isinstance(canonical_smiles, str) or len(canonical_smiles) > 2000:
            raise ValueError("canonical_smiles 必须为最多 2000 个字符的字符串。")
        with closing(self._connect()) as connection:
            rows = connection.execute("""SELECT o.dataset,o.source_id,o.target,o.value,o.unit,o.quality,d.source_url
                FROM observations AS o JOIN molecules AS m ON o.molecule_id=m.id
                JOIN datasets AS d ON o.dataset=d.id WHERE m.canonical_smiles=?
                ORDER BY o.target,o.dataset,o.source_id""", (canonical_smiles,))
            return [dict(row) for row in rows]
