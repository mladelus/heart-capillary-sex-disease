"""Minimal reader for R .rds files (XDR, version 2/3). Returns dicts / pandas DataFrames / numpy arrays."""
import gzip, bz2, lzma, struct
import numpy as np, pandas as pd

NA_INT = -2147483648


class _R:
    def __init__(self, b):
        self.b = b; self.i = 0; self.refs = []

    def int(self):
        v = struct.unpack_from(">i", self.b, self.i)[0]; self.i += 4; return v

    def length(self):
        n = self.int()
        if n == -1:
            hi = self.int(); lo = self.int(); n = (hi << 32) + (lo & 0xFFFFFFFF)
        return n

    def item(self):
        flags = self.int()
        t = flags & 0xFF; has_attr = bool(flags & 0x200); has_tag = bool(flags & 0x400)
        if t in (254, 253, 252, 251, 250, 249, 248):
            return None
        if t == 255:
            idx = flags >> 8
            if idx == 0: idx = self.int()
            return self.refs[idx - 1]
        if t == 1:  # symbol
            s = self.item(); self.refs.append(s); return s
        if t in (247, 246, 245):  # namespace / package / persist
            self.int(); n = self.int(); s = [self.item() for _ in range(n)]
            self.refs.append(("ns", s)); return ("ns", s)
        if t == 4:  # environment
            env = {"__env__": True}; self.refs.append(env)
            self.int(); self.item(); self.item(); self.item(); self.item(); return env
        if t in (2, 3, 5, 6, 17, 240, 243):  # pairlist-like
            out = []
            while True:
                attr = self.item() if has_attr or t in (240, 243) and False else (self.item() if has_attr else None)
                tag = self.item() if has_tag else None
                car = self.item()
                out.append((tag, car))
                nflags_pos = self.i
                nflags = self.int()
                nt = nflags & 0xFF
                if nt == 254:
                    break
                if nt != 2 or t != 2:
                    self.i = nflags_pos; out.append((None, self.item())); break
                has_attr = bool(nflags & 0x200); has_tag = bool(nflags & 0x400)
            return out
        if t == 9:  # CHARSXP
            n = self.int()
            if n == -1: return None
            s = self.b[self.i:self.i + n]; self.i += n
            return s.decode("utf-8", "replace")
        if t == 238:  # ALTREP
            info = self.item(); state = self.item(); attr = self.item()
            cls = info[0][1]
            if cls in ("compact_intseq", "compact_realseq"):
                n, start, step = state["v"] if isinstance(state, dict) else state
                v = start + step * np.arange(int(n))
                val = v.astype(int) if cls == "compact_intseq" else v.astype(float)
            elif cls == "deferred_string":
                arg = state[0][1]; arg = arg["v"] if isinstance(arg, dict) else arg
                val = np.array([None if (isinstance(x, float) and np.isnan(x)) else (str(int(x)) if float(x).is_integer() else repr(float(x))) for x in np.asarray(arg, dtype=float)], dtype=object)
            elif cls.startswith("wrap_"):
                val = state[0] if isinstance(state, list) else state
                if isinstance(val, dict) and "v" in val and "a" in val and len(val) == 2: return self._wrap(val["v"], {**val["a"], **self._attrs(attr)})
            else:
                raise ValueError("altrep " + cls)
            return self._wrap(val, self._attrs(attr))
        if t == 10:
            n = self.length(); v = np.frombuffer(self.b, ">i4", n, self.i).astype(float); self.i += 4 * n
            v[v == NA_INT] = np.nan; val = v
        elif t == 13:
            n = self.length(); val = np.frombuffer(self.b, ">i4", n, self.i).astype(np.int64); self.i += 4 * n
        elif t == 14:
            n = self.length(); val = np.frombuffer(self.b, ">f8", n, self.i).astype(float); self.i += 8 * n
        elif t == 15:
            n = self.length(); val = np.frombuffer(self.b, ">f8", 2 * n, self.i).astype(float); self.i += 16 * n
        elif t == 16:
            n = self.length(); val = np.array([self.item() for _ in range(n)], dtype=object)
        elif t in (19, 20):
            n = self.length(); val = [self.item() for _ in range(n)]
        elif t == 24:
            n = self.length(); val = self.b[self.i:self.i + n]; self.i += n
        elif t == 22:  # external pointer
            self.refs.append(None); self.item(); self.item(); val = None
        elif t == 25:  # S4 object: keep its slots
            attrs = self._attrs(self.item()) if has_attr else {}
            return {"__S4__": True, **attrs}
        else:
            raise ValueError(f"type {t} at {self.i}")
        attrs = self._attrs(self.item()) if has_attr else {}
        return self._wrap(val, attrs)

    @staticmethod
    def _attrs(pl):
        if not pl: return {}
        return {k: v for k, v in pl if k is not None}

    @staticmethod
    def _wrap(val, attrs):
        if not attrs: return val
        cls = attrs.get("class"); cls = list(cls) if cls is not None and not isinstance(cls, dict) else []
        if isinstance(val, np.ndarray) and "levels" in attrs and val.dtype != object:
            lev = list(attrs["levels"])
            return pd.Categorical.from_codes(np.where(val == NA_INT, -1, val - 1).astype(int), categories=lev, ordered="ordered" in cls)
        if isinstance(val, list):
            names = attrs.get("names")
            if "data.frame" in cls:
                cols = {}
                for k, v in zip(names, val):
                    if isinstance(v, np.ndarray) and v.dtype == np.int64:
                        v = np.where(v == NA_INT, np.nan, v) if (v == NA_INT).any() else v
                    cols[k] = v if not isinstance(v, (list, dict)) else pd.Series(list(v) if isinstance(v, list) else [v], dtype=object)
                try:
                    return pd.DataFrame(cols)
                except Exception:
                    return {"__df_failed__": True, **{k: v for k, v in zip(names, val)}}
            if names is not None:
                d = {}
                for j, (k, v) in enumerate(zip(names, val)):
                    d[k if k else f"_{j}"] = v
                if any(a not in ("names", "class") for a in attrs): d["__attrs__"] = {a: b for a, b in attrs.items() if a != "names"}
                return d
            return val
        if isinstance(val, np.ndarray):
            if "dim" in attrs and val.dtype != object:
                dim = [int(x) for x in attrs["dim"]]
                m = val.reshape(dim, order="F")
                dn = attrs.get("dimnames")
                if len(dim) == 2 and dn is not None:
                    return pd.DataFrame(m, index=list(dn[0]) if dn[0] is not None else None, columns=list(dn[1]) if dn[1] is not None else None)
                return m
            if "names" in attrs and val.dtype != object:
                return pd.Series(val, index=list(attrs["names"]))
            if "names" in attrs:
                return pd.Series(list(val), index=list(attrs["names"]), dtype=object)
            if len(attrs) and set(attrs) - {"names"}:
                return {"v": val, "a": attrs} if False else val
        return val


def read_rds(path):
    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b": raw = gzip.decompress(raw)
    elif raw[:3] == b"BZh": raw = bz2.decompress(raw)
    elif raw[:6] == b"\xfd7zXZ\x00": raw = lzma.decompress(raw)
    assert raw[:2] == b"X\n", raw[:10]
    r = _R(raw); r.i = 2
    ver = r.int(); r.int(); r.int()
    if ver >= 3:
        n = r.int(); r.i += n
    return r.item()


def show(x, ind=0, name=""):
    pad = " " * ind
    if isinstance(x, pd.DataFrame): print(f"{pad}{name}: DataFrame {x.shape} {list(x.columns)[:40]}")
    elif isinstance(x, dict):
        print(f"{pad}{name}: dict")
        for k, v in x.items():
            if k != "__attrs__": show(v, ind + 2, k)
    elif isinstance(x, list): print(f"{pad}{name}: list[{len(x)}]"); [show(v, ind + 2, str(j)) for j, v in enumerate(x[:2])]
    elif isinstance(x, (np.ndarray, pd.Series, pd.Categorical)): print(f"{pad}{name}: {type(x).__name__}[{len(x)}]")
    else: print(f"{pad}{name}: {type(x).__name__} {str(x)[:60]}")
