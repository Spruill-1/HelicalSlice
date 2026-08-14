"""Inject the thumbnail set Bambu Connect requires into headless-sliced .gcode.3mf files.

Renders honest top-view images from Metadata/plate_1.json object bboxes:
plate_1.png (512), plate_1_small.png (128), plate_no_light_1.png, top_1.png, pick_1.png,
and registers them in the model_settings.config plate block.
"""
import zipfile, json, os, re
from PIL import Image, ImageDraw

SLICED = r"C:\Users\david\source\HelicalSlice\test_models\tpu_suite\sliced"
BED_W, BED_H = 350.0, 320.0
SZ = 512

FILES = ["PlateA_T1-T4", "T5_seam_scarf0", "T5_seam_scarf6", "T5_seam_scarf12", "T5_trio_fmb"]
TUBE_COLORS = [(46, 134, 222), (235, 77, 75), (38, 222, 129), (247, 159, 31), (162, 155, 254)]

def bed_to_px(x, y):
    # bed origin bottom-left -> image top-left; keep aspect via uniform scale
    s = SZ / max(BED_W, BED_H)
    ox = (SZ - BED_W * s) / 2
    oy = (SZ - BED_H * s) / 2
    return ox + x * s, SZ - (oy + y * s)

def render(objs):
    plate = Image.new("RGB", (SZ, SZ), (38, 42, 48))
    d = ImageDraw.Draw(plate)
    # bed outline + grid
    x0, y0 = bed_to_px(0, 0); x1, y1 = bed_to_px(BED_W, BED_H)
    d.rectangle([x1 if x1 < x0 else x0, y1 if y1 < y0 else y0, max(x0, x1), max(y0, y1)], fill=(52, 58, 66), outline=(90, 98, 108), width=2)
    for gx in range(0, int(BED_W) + 1, 50):
        a = bed_to_px(gx, 0); b = bed_to_px(gx, BED_H); d.line([a, b], fill=(64, 70, 78), width=1)
    for gy in range(0, int(BED_H) + 1, 50):
        a = bed_to_px(0, gy); b = bed_to_px(BED_W, gy); d.line([a, b], fill=(64, 70, 78), width=1)

    top  = Image.new("RGBA", (SZ, SZ), (0, 0, 0, 0));  dt = ImageDraw.Draw(top)
    pick = Image.new("RGBA", (SZ, SZ), (0, 0, 0, 0));  dp = ImageDraw.Draw(pick)

    for i, o in enumerate(objs):
        bx0, by0, bx1, by1 = o["bbox"]
        col = TUBE_COLORS[i % len(TUBE_COLORS)]
        pa = bed_to_px(bx0, by0); pb = bed_to_px(bx1, by1)
        box = [min(pa[0], pb[0]), min(pa[1], pb[1]), max(pa[0], pb[0]), max(pa[1], pb[1])]
        wall = max(3, int((box[2] - box[0]) * 0.12))
        for img_d, c in ((d, col), (dt, col + (255,))):
            img_d.ellipse(box, outline=c, width=wall)
        # pick: object index encoded in the red channel (i+1), flat filled silhouette
        dp.ellipse(box, fill=(i + 1, 0, 0, 255))
    return plate, top, pick

def png_bytes(img):
    from io import BytesIO
    b = BytesIO(); img.save(b, "PNG"); return b.getvalue()

for name in FILES:
    path = os.path.join(SLICED, name + ".gcode.3mf")
    zin = zipfile.ZipFile(path)
    pj = json.loads(zin.read("Metadata/plate_1.json"))
    plate, top, pick = render(pj["bbox_objects"])
    small = plate.resize((128, 128), Image.LANCZOS)

    ms = zin.read("Metadata/model_settings.config").decode("utf-8")
    if "thumbnail_file" not in ms:
        ms = ms.replace('<metadata key="gcode_file" value="Metadata/plate_1.gcode"/>',
                        '<metadata key="gcode_file" value="Metadata/plate_1.gcode"/>\n'
                        '    <metadata key="thumbnail_file" value="Metadata/plate_1.png"/>\n'
                        '    <metadata key="thumbnail_no_light_file" value="Metadata/plate_no_light_1.png"/>\n'
                        '    <metadata key="top_file" value="Metadata/top_1.png"/>\n'
                        '    <metadata key="pick_file" value="Metadata/pick_1.png"/>', 1)

    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            if info.filename == "Metadata/model_settings.config":
                zout.writestr(info.filename, ms)
            else:
                zout.writestr(info, zin.read(info.filename))
        zout.writestr("Metadata/plate_1.png", png_bytes(plate))
        zout.writestr("Metadata/plate_1_small.png", png_bytes(small))
        zout.writestr("Metadata/plate_no_light_1.png", png_bytes(plate))
        zout.writestr("Metadata/top_1.png", png_bytes(top))
        zout.writestr("Metadata/pick_1.png", png_bytes(pick))
    zin.close()
    os.replace(tmp, path)
    z = zipfile.ZipFile(path)
    pngs = [n for n in z.namelist() if n.endswith(".png")]
    print(f"{name:22s} injected: {pngs}")
