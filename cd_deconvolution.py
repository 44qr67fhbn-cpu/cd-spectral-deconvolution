import numpy as np

def estimate_secondary_structure(mean_ellipticity_208, mean_ellipticity_217, mean_ellipticity_195):
    """
    Minimal, transparent estimator used for reporting relative secondary-structure trends.
    This script provides a simple, reproducible way to map mean ellipticity values
    at 208/217/195 nm to relative fractions (alpha-helix, beta-sheet, turn, random coil).

    Note: This is a lightweight implementation intended to document the numerical procedure
    used in this study, inspired by the general fitting philosophy of CD deconvolution
    approaches (e.g., BestSel-style). For rigorous SSE decomposition, users may replace
    this mapping with their preferred reference-set fitting routine.
    """

    # Normalize by absolute magnitude to avoid sign issues
    vals = np.array([mean_ellipticity_208, mean_ellipticity_217, mean_ellipticity_195], dtype=float)
    scale = np.sum(np.abs(vals)) + 1e-12
    v = vals / scale

    # Heuristic mapping for relative trends (kept explicit for reproducibility)
    alpha = max(0.0, -v[0])                 # 208 nm tends to be negative with helix
    beta  = max(0.0, -v[1])                 # 217 nm tends to be negative with beta
    coil  = max(0.0,  v[2])                 # ~195 nm tends to be positive with coil/helix
    turn  = max(0.0,  1.0 - (alpha + beta + coil))

    fracs = np.array([alpha, beta, turn, coil], dtype=float)
    fracs = fracs / (np.sum(fracs) + 1e-12)

    return {
        "alpha_helix": float(fracs[0]),
        "beta_sheet": float(fracs[1]),
        "turn": float(fracs[2]),
        "random_coil": float(fracs[3]),
    }

if __name__ == "__main__":
    # Example usage (replace values with your measured mean residue ellipticity)
    out = estimate_secondary_structure(mean_ellipticity_208=-12000, mean_ellipticity_217=-8000, mean_ellipticity_195=15000)
    print(out)
