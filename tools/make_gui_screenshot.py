"""Capture simulator/cassegrain_gui.py into
docs/figures/cassegrain_gui_screenshot.png (Windows)."""

import ctypes
ctypes.windll.shcore.SetProcessDpiAwareness(1)   # physical-pixel coordinates

import sys
import tkinter as tk
from pathlib import Path
from PIL import ImageGrab

ROOT = Path(__file__).resolve().parent.parent
OUT_PNG = ROOT / "docs" / "figures" / "cassegrain_gui_screenshot.png"
sys.path.insert(0, str(ROOT / "simulator"))
import cassegrain_gui as G


def main():
    root = tk.Tk()
    root.geometry("+40+40")
    app = G.App(root)
    # representative state: polychromatic donut of a small extended source
    app.sliders["defocus"].set(16.0)
    for key, value in (("lower_um", "0.5225"), ("upper_um", "0.5775")):
        app.entries[key].delete(0, tk.END)
        app.entries[key].insert(0, value)
    app._apply_params()
    app.sliders["source"].set(2.0)
    root.update_idletasks()
    root.update()
    app._update_image()                          # synchronous final render
    root.attributes("-topmost", True)
    root.update_idletasks()
    root.update()

    x, y = root.winfo_rootx(), root.winfo_rooty()
    w, h = root.winfo_width(), root.winfo_height()
    img = ImageGrab.grab(bbox=(x, y, x + w, y + h))
    img.save(OUT_PNG)
    print("saved cassegrain_gui_screenshot.png", img.size)
    root.destroy()


if __name__ == "__main__":
    main()
