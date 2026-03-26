---
name: Do not read large dachuang images
description: Never read marked_*.png, raw_images/*, or any generated marked_ images in dachuang directory
type: feedback
---

Do not read these files in `/home/gyqyan/zexuan/ecg1d_cpp/dachuang/`:
- `marked_0ugL.png`, `marked_5ugL.png`, `marked_20ugL.png`
- Anything in `raw_images/`
- Any newly generated `marked_*` images

**Why:** User explicitly forbids it — these are large images (12MB, 3536x3536) that waste context. When code generates `marked_` images, just inform the user they're ready and let them view manually.

**How to apply:** When working in the dachuang directory, never use the Read tool on these files. Only read images the user explicitly asks to see (e.g., `_subregion_*.png` analysis charts). If analysis code produces `marked_` output, report completion without reading the file.
