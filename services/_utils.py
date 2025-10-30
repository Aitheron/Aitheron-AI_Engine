import pandas as pd

def drop_cols(
    df: pd.DataFrame,
    keep: list[str] = None,
    also_drop: list[str] = None,
    prefix: str | list[str] | tuple[str, ...] = "is",
    inplace: bool = False,
) -> pd.DataFrame:
    keep = set(keep or [])
    cols = df.columns

    if isinstance(prefix, (list, tuple, set)):
        pref = tuple(str(p).lower() for p in prefix)
        mask = cols.str.lower().str.startswith(pref)
    else:
        p = str(prefix).lower()
        mask = cols.str.lower().str.startswith(p)

    to_drop = [c for c in cols[mask] if c not in keep]

    if also_drop:
        for c in also_drop:
            if c not in keep:
                to_drop.append(c)

    to_drop = list(dict.fromkeys(to_drop))  # dedup
    if inplace:
        df.drop(columns=to_drop, errors="ignore", inplace=True)
        return df
    return df.drop(columns=to_drop, errors="ignore")

def reorder_columns(df, moves, after=True):
    # moves: lista de pares (col, where_col)
    cols = list(df.columns)
    for col, where in moves:
        if col in cols and where in cols and col != where:
            cols.remove(col)
            idx = cols.index(where)
            cols.insert(idx + (1 if after else 0), col)
    return df.loc[:, cols]
