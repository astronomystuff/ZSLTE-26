import requests
import numpy as np

# Horizons API endpoint
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

def fetch_state(target, epoch="2451545.0"):
    """
    Fetch barycentric state vector from JPL Horizons.
    Returns dict with x,y,z,vx,vy,vz in AU and AU/day.
    """

    params = {
        "format": "json",
        "COMMAND": target,
        "CENTER": "500@0",          # Solar System Barycenter
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "OUT_UNITS": "AU-D",
        "REF_PLANE": "ECLIPTIC",
        "REF_SYSTEM": "J2000",
        "TLIST": epoch
    }

    r = requests.get(HORIZONS_URL, params=params)
    data = r.json()

    try:
        vec = data["result"].split("$$SOE")[1].split("$$EOE")[0].strip().split("\n")[1]
    except Exception:
        raise RuntimeError(f"Failed to parse Horizons response for {target}")

    fields = vec.split()

    # Horizons returns:
    # x y z vx vy vz (AU and AU/day)
    x, y, z = map(float, fields[2:5])
    vx, vy, vz = map(float, fields[5:8])

    return {
        "x": x,
        "y": y,
        "z": z,
        "vx": vx,
        "vy": vy,
        "vz": vz
    }
