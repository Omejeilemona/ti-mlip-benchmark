"""Minimal ray-cast sphere renderer for atomistic snapshots (NumPy only, no GPU)."""
import numpy as np

def rotation(view):
    d = np.asarray(view, float); d /= np.linalg.norm(d)
    up = np.array([0, 0, 1.0]) if abs(d[2]) < 0.9 else np.array([0, 1.0, 0])
    x = np.cross(up, -d); x /= np.linalg.norm(x); y = np.cross(-d, x)
    return np.array([x, y, -d])

def render(pos, colors, radius=1.3, fov=40.0, size=600, view=(0.35, -0.55, -0.76), light=(-0.5, 0.6, 1.0),
           background=(1, 1, 1)):
    R = rotation(view); p = pos @ R.T
    img = np.ones((size, size, 3)) * np.asarray(background, float)
    zbuf = np.full((size, size), -np.inf)
    scale = size / (2 * fov)
    L = np.asarray(light, float); L /= np.linalg.norm(L); H = L + np.array([0, 0, 1.0]); H /= np.linalg.norm(H)
    zmin, zmax = p[:, 2].min(), p[:, 2].max()
    rpx = radius * scale
    for i in np.argsort(p[:, 2]):
        cx, cy = (p[i, 0] + fov) * scale, (fov - p[i, 1]) * scale
        x0, x1 = int(max(cx - rpx - 1, 0)), int(min(cx + rpx + 2, size))
        y0, y1 = int(max(cy - rpx - 1, 0)), int(min(cy + rpx + 2, size))
        if x0 >= x1 or y0 >= y1: continue
        yy, xx = np.mgrid[y0:y1, x0:x1]
        dx, dy = (xx + 0.5 - cx) / rpx, (cy - yy - 0.5) / rpx
        r2 = dx * dx + dy * dy; inside = r2 < 1.0
        if not inside.any(): continue
        nz = np.sqrt(np.clip(1 - r2, 0, 1)); z = p[i, 2] + nz * radius
        upd = inside & (z > zbuf[y0:y1, x0:x1])
        if not upd.any(): continue
        n = np.stack([dx, dy, nz], -1)
        lam = np.clip(n @ L, 0, 1); spec = np.clip(n @ H, 0, 1) ** 40
        depth = 0.75 + 0.25 * (p[i, 2] - zmin) / max(zmax - zmin, 1e-9)
        rim = np.clip((1 - r2) / 0.12, 0, 1) ** 0.5
        c = np.asarray(colors[i], float)
        shade = (0.30 + 0.70 * lam)[..., None] * c * depth * rim[..., None] + 0.45 * spec[..., None] * rim[..., None]
        img[y0:y1, x0:x1][upd] = np.clip(shade[upd], 0, 1)
        zbuf[y0:y1, x0:x1][upd] = z[upd]
    return img
