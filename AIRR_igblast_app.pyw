# -*- coding: utf-8 -*-
import csv
import math
import os
import shutil
import subprocess
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import ctypes
from xml.sax.saxutils import escape
import tkinter as tk
from tkinter import filedialog, messagebox

APP_DIR = Path(__file__).resolve().parent
RESULT_DIR = APP_DIR / "result_AIRR_outfmat"
WORK_ROOT = Path(os.environ.get("LOCALAPPDATA", str(APP_DIR))) / "RGMergedFastaIgblastAirrTsv"
# Full path to reference data (can include non-ASCII).
REF_DIR_FULL = Path(r"C:\Users\Yohei Funakoshi\Desktop\IgBlast用参照データ")
IGBLAST = Path(r"C:\Program Files\NCBI\igblast-1.21.0\bin\igblastn.exe")
REF_DIR_USE = None
NUM_THREADS = max(1, min(4, os.cpu_count() or 1))
FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

KERNEL32 = ctypes.windll.kernel32
KERNEL32.GetShortPathNameW.argtypes = [ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint]
KERNEL32.GetShortPathNameW.restype = ctypes.c_uint


def short_path(path_str):
    size = max(260, len(path_str) + 1)
    while size <= 32768:
        buf = ctypes.create_unicode_buffer(size)
        ret = KERNEL32.GetShortPathNameW(path_str, buf, size)
        if ret == 0:
            return path_str
        if ret < size:
            return buf.value
        size = ret + 1
    return path_str


def get_ref_dir():
    global REF_DIR_USE
    if REF_DIR_USE is not None:
        return REF_DIR_USE

    ref_short = Path(short_path(str(REF_DIR_FULL)))
    if (ref_short / "db").exists():
        REF_DIR_USE = ref_short
        return REF_DIR_USE

    # Create an ASCII junction under LocalAppData if 8.3 short path isn't available.
    junction = WORK_ROOT / "refdata"
    if not junction.exists():
        try:
            junction.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                ["cmd", "/c", "mklink", "/J", str(junction), str(REF_DIR_FULL)],
                capture_output=True,
                text=True,
            )
        except Exception:
            pass
    if junction.exists():
        REF_DIR_USE = junction
        return REF_DIR_USE

    REF_DIR_USE = REF_DIR_FULL
    return REF_DIR_USE


def collect_missing_prereqs():
    missing = []
    if not IGBLAST.exists():
        missing.append(str(IGBLAST))
    ref_dir = get_ref_dir()
    db_dir = ref_dir / "db"
    aux_file = ref_dir / "optional_file" / "human_gl.aux"
    internal_file = ref_dir / "internal_data" / "human" / "human.ndm.imgt"
    for name in ("IMGT_IGHV.imgt.nsq", "IMGT_IGHD.imgt.nsq", "IMGT_IGHJ.imgt.nsq"):
        if not (db_dir / name).exists():
            missing.append(str(db_dir / name))
    if not aux_file.exists():
        missing.append(str(aux_file))
    if not internal_file.exists():
        missing.append(str(internal_file))
    return missing


def check_prereq():
    missing = collect_missing_prereqs()
    if missing:
        msg = "Missing files:\n" + "\n".join(missing)
        messagebox.showerror("Missing files", msg)
        return False
    return True


def inspect_fasta(in_path):
    if in_path.suffix.lower() not in FASTA_EXTENSIONS:
        raise ValueError("Input file is not .fasta, .fa, or .fna.")

    records = 0
    bases = 0
    first_id = ""
    allowed = set("ACGTRYSWKMBDHVN.-")
    with in_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records += 1
                if not first_id:
                    first_id = line[1:].strip()[:120]
                continue
            if records == 0:
                raise ValueError(f"Line {line_no}: first non-empty line is not a FASTA header.")
            seq = "".join(line.split()).upper()
            invalid = sorted(set(seq) - allowed)
            if invalid:
                raise ValueError(f"Line {line_no}: invalid FASTA character(s): {', '.join(invalid)}")
            bases += len(seq)

    if records == 0:
        raise ValueError("No FASTA records found.")
    if bases == 0:
        raise ValueError("No sequence bases found.")
    return {"records": records, "bases": bases, "first_id": first_id}


def make_work_dir():
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    work_dir = WORK_ROOT / "work" / stamp
    work_dir.mkdir(parents=True, exist_ok=True)
    return work_dir


def output_base_name(in_path):
    stem = in_path.stem
    base = stem
    for marker in (
        "_trimThenMerge_",
        "_mergeThenTrim_",
        "_withoutTrim_",
        "_presto_assemble-pass",
        "_assemble-pass",
        "_extendedFrags",
    ):
        if marker in base:
            base = base.split(marker, 1)[0]
            break

    labels = []
    lower_stem = stem.lower()
    if "filtered" in lower_stem or "no_suspect" in lower_stem:
        labels.append("filtered")
    if not base:
        base = stem[:60]
    if labels:
        base = base + "." + ".".join(labels)
    return base[:90]


def is_truthy(value):
    return str(value or "").strip().lower() in {"t", "true", "1", "yes", "y"}


def gene_set(call):
    genes = []
    for part in str(call or "").replace(";", ",").split(","):
        gene = part.strip()
        if not gene:
            continue
        gene = gene.split()[0]
        gene = gene.split("*", 1)[0]
        if gene:
            genes.append(gene)
    return ",".join(sorted(set(genes)))


def canonical_junction_aa(seq):
    seq = str(seq or "").strip().upper()
    if not seq:
        return False, "missing_junction_aa"
    if "*" in seq:
        return False, "stop_codon"
    if len(seq) < 5 or len(seq) > 40:
        return False, "junction_aa_length"
    if not seq.startswith("C"):
        return False, "junction_aa_no_C_start"
    if seq[-1] not in {"W", "F"}:
        return False, "junction_aa_no_WF_end"
    if any(ch not in AA_ALPHABET for ch in seq):
        return False, "invalid_junction_aa"
    return True, ""


def safe_float(value):
    try:
        if value is None or value == "":
            return None
        return float(value)
    except ValueError:
        return None


def summarize_airr(airr_path, counts_tsv_path, xlsx_path):
    total_rows = 0
    included_rows = 0
    productive_true_rows = 0
    excluded = Counter()
    groups = defaultdict(
        lambda: {
            "read_count": 0,
            "productive_true_count": 0,
            "v_identity_values": [],
            "representative_sequence_id": "",
            "representative_v_call": "",
            "representative_d_call": "",
            "representative_j_call": "",
            "representative_locus": "",
        }
    )

    with airr_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("AIRR TSV has no header.")

        required = {"sequence_id", "v_call", "j_call", "junction_aa"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError("AIRR TSV is missing required columns: " + ", ".join(missing))

        for row in reader:
            total_rows += 1
            locus = (row.get("locus") or "").strip()
            if locus and locus != "IGH":
                excluded["non_IGH"] += 1
                continue

            v_set = gene_set(row.get("v_call"))
            j_set = gene_set(row.get("j_call"))
            if not v_set:
                excluded["missing_v_call"] += 1
                continue
            if not j_set:
                excluded["missing_j_call"] += 1
                continue
            if not is_truthy(row.get("productive")):
                excluded["not_productive"] += 1
                continue

            junction_aa = (row.get("junction_aa") or "").strip().upper()
            ok, reason = canonical_junction_aa(junction_aa)
            if not ok:
                excluded[reason] += 1
                continue

            key = (v_set, j_set, junction_aa)
            group = groups[key]
            group["read_count"] += 1
            included_rows += 1
            if is_truthy(row.get("productive")):
                group["productive_true_count"] += 1
                productive_true_rows += 1
            v_identity = safe_float(row.get("v_identity"))
            if v_identity is not None:
                group["v_identity_values"].append(v_identity)
            if not group["representative_sequence_id"]:
                group["representative_sequence_id"] = row.get("sequence_id", "")
                group["representative_v_call"] = row.get("v_call", "")
                group["representative_d_call"] = row.get("d_call", "")
                group["representative_j_call"] = row.get("j_call", "")
                group["representative_locus"] = locus

    headers = [
        "unique_v_gene_set",
        "unique_j_gene_set",
        "junction_aa",
        "junction_aa_length",
        "read_count",
        "productive_true_count",
        "productive_true_rate",
        "avg_v_identity",
        "min_v_identity",
        "max_v_identity",
        "representative_sequence_id",
        "representative_v_call",
        "representative_d_call",
        "representative_j_call",
        "locus",
    ]

    count_rows = []
    for (v_set, j_set, junction_aa), group in groups.items():
        values = group["v_identity_values"]
        avg_v = sum(values) / len(values) if values else ""
        min_v = min(values) if values else ""
        max_v = max(values) if values else ""
        read_count = group["read_count"]
        prod_count = group["productive_true_count"]
        count_rows.append(
            [
                v_set,
                j_set,
                junction_aa,
                len(junction_aa),
                read_count,
                prod_count,
                prod_count / read_count if read_count else 0,
                avg_v,
                min_v,
                max_v,
                group["representative_sequence_id"],
                group["representative_v_call"],
                group["representative_d_call"],
                group["representative_j_call"],
                group["representative_locus"],
            ]
        )
    count_rows.sort(key=lambda row: (-row[4], row[0], row[1], row[2]))

    with counts_tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(count_rows)

    summary_rows = [
        ["Metric", "Value"],
        ["AIRR TSV", str(airr_path)],
        ["Generated", datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        ["Total AIRR rows", total_rows],
        ["Rows included in counts", included_rows],
        ["Unique clonotypes", len(count_rows)],
        ["Productive true rows in counts", productive_true_rows],
        ["Counts rule", "IGH + productive + V call + J call + canonical junction_aa"],
    ]
    excluded_rows = [["Exclude reason", "Rows"]] + [[reason, count] for reason, count in excluded.most_common()]
    write_xlsx(
        xlsx_path,
        [
            ("Summary", summary_rows),
            ("Counts", [headers] + count_rows),
            ("Excluded", excluded_rows),
        ],
    )

    return {
        "total_rows": total_rows,
        "included_rows": included_rows,
        "unique_clonotypes": len(count_rows),
        "productive_true_rows": productive_true_rows,
        "excluded": excluded,
        "counts_tsv": counts_tsv_path,
        "xlsx": xlsx_path,
    }


def xml_text(value):
    text = "" if value is None else str(value)
    text = "".join(ch for ch in text if ch == "\t" or ch == "\n" or ch == "\r" or ord(ch) >= 32)
    return escape(text, {"'": "&apos;", '"': "&quot;"})


def column_name(index):
    name = ""
    while index:
        index, rem = divmod(index - 1, 26)
        name = chr(65 + rem) + name
    return name


def cell_xml(row_idx, col_idx, value, header=False):
    ref = f"{column_name(col_idx)}{row_idx}"
    style = ' s="1"' if header else ""
    if value is None or value == "":
        return f'<c r="{ref}"{style}/>'
    if isinstance(value, bool):
        return f'<c r="{ref}" t="b"{style}><v>{1 if value else 0}</v></c>'
    if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        return f'<c r="{ref}"{style}><v>{value}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{style}><is><t>{xml_text(value)}</t></is></c>'


def sheet_xml(rows, freeze_header=True):
    max_cols = max((len(row) for row in rows), default=1)
    max_rows = max(len(rows), 1)
    col_widths = []
    for col_idx in range(max_cols):
        max_len = 8
        for row in rows[:1000]:
            if col_idx < len(row):
                max_len = max(max_len, len(str(row[col_idx])))
        col_widths.append(min(max(max_len + 2, 10), 45))

    cols_xml = "".join(
        f'<col min="{i}" max="{i}" width="{width}" customWidth="1"/>'
        for i, width in enumerate(col_widths, start=1)
    )
    pane_xml = ""
    if freeze_header and rows:
        pane_xml = (
            '<sheetViews><sheetView workbookViewId="0">'
            '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
            '</sheetView></sheetViews>'
        )
    row_xml = []
    for row_idx, row in enumerate(rows, start=1):
        cells = "".join(cell_xml(row_idx, col_idx, value, header=(row_idx == 1)) for col_idx, value in enumerate(row, start=1))
        row_xml.append(f'<row r="{row_idx}">{cells}</row>')
    auto_filter = f'<autoFilter ref="A1:{column_name(max_cols)}{max_rows}"/>' if max_rows > 1 else ""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"{pane_xml}<cols>{cols_xml}</cols><sheetData>{''.join(row_xml)}</sheetData>{auto_filter}</worksheet>"
    )


def safe_sheet_name(name, used):
    invalid = set(r'[]:*?/\\')
    cleaned = "".join("_" if ch in invalid else ch for ch in name).strip() or "Sheet"
    cleaned = cleaned[:31]
    base = cleaned
    suffix = 1
    while cleaned in used:
        suffix_text = f"_{suffix}"
        cleaned = (base[: 31 - len(suffix_text)] + suffix_text)[:31]
        suffix += 1
    used.add(cleaned)
    return cleaned


def write_xlsx(path, sheets):
    used = set()
    named_sheets = [(safe_sheet_name(name, used), rows) for name, rows in sheets]
    sheet_overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for i in range(1, len(named_sheets) + 1)
    )
    workbook_sheets = "".join(
        f'<sheet name="{xml_text(name)}" sheetId="{i}" r:id="rId{i}"/>'
        for i, (name, _rows) in enumerate(named_sheets, start=1)
    )
    workbook_rels = "".join(
        f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
        for i in range(1, len(named_sheets) + 1)
    )
    styles_rid = len(named_sheets) + 1
    now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
            '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
            f"{sheet_overrides}</Types>",
        )
        zf.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
            '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>'
            '</Relationships>',
        )
        zf.writestr(
            "docProps/core.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" '
            'xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '<dc:creator>RG merged FASTA AIRR app</dc:creator>'
            f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
            f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>'
            '</cp:coreProperties>',
        )
        zf.writestr(
            "docProps/app.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
            'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
            '<Application>Microsoft Excel</Application></Properties>',
        )
        zf.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f"<sheets>{workbook_sheets}</sheets></workbook>",
        )
        zf.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{workbook_rels}<Relationship Id="rId{styles_rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
            '</Relationships>',
        )
        zf.writestr(
            "xl/styles.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font></fonts>'
            '<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill></fills>'
            '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/></cellXfs>'
            '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
            '</styleSheet>',
        )
        for i, (_name, rows) in enumerate(named_sheets, start=1):
            zf.writestr(f"xl/worksheets/sheet{i}.xml", sheet_xml(rows))


def create_analysis_outputs(airr_path):
    suffix = ".igblast.airr.tsv"
    if airr_path.name.endswith(suffix):
        base_name = airr_path.name[: -len(suffix)]
    else:
        base_name = airr_path.stem
    counts_tsv_path = airr_path.parent / f"{base_name}.airr_counts.tsv"
    xlsx_path = airr_path.parent / f"{base_name}.airr_counts.xlsx"
    return summarize_airr(airr_path, counts_tsv_path, xlsx_path)


def run_igblast(input_path, output_dir, log):
    if not check_prereq():
        return
    if not input_path:
        messagebox.showwarning("Input", "Please select a FASTA file.")
        return
    in_path = Path(input_path)
    if not in_path.exists():
        messagebox.showerror("Input", "Input file not found.")
        return
    try:
        fasta_info = inspect_fasta(in_path)
    except Exception as exc:
        messagebox.showerror("Input FASTA", str(exc))
        return

    out_dir = Path(output_dir) if output_dir else RESULT_DIR
    out_dir.mkdir(exist_ok=True)
    out_name = output_base_name(in_path) + ".igblast.airr.tsv"
    out_path = out_dir / out_name
    work_dir = make_work_dir()
    work_query = work_dir / "query.fasta"
    work_out = work_dir / "igblast.airr.tsv"
    shutil.copy2(in_path, work_query)

    env = os.environ.copy()
    ref_dir = get_ref_dir()
    ref_short = short_path(str(ref_dir))
    db_dir = ref_dir / "db"
    aux_file = ref_dir / "optional_file" / "human_gl.aux"
    env["IGDATA"] = ref_short

    args = [
        str(IGBLAST),
        "-query", short_path(str(work_query)),
        "-germline_db_V", short_path(str(db_dir / "IMGT_IGHV.imgt")),
        "-germline_db_D", short_path(str(db_dir / "IMGT_IGHD.imgt")),
        "-germline_db_J", short_path(str(db_dir / "IMGT_IGHJ.imgt")),
        "-auxiliary_data", short_path(str(aux_file)),
        "-domain_system", "imgt",
        "-organism", "human",
        "-ig_seqtype", "Ig",
        "-outfmt", "19",
        "-num_threads", str(NUM_THREADS),
        "-out", short_path(str(work_out)),
    ]

    log(f"Input FASTA: {fasta_info['records']:,} records / {fasta_info['bases']:,} bases")
    log(f"Work folder: {work_dir}")
    log(f"Running IgBLAST with {NUM_THREADS} thread(s)...")
    try:
        proc = subprocess.run(args, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
    except Exception as exc:
        messagebox.showerror("Error", f"Failed to run: {exc}")
        return

    if proc.stdout:
        log(proc.stdout.strip())
    if proc.stderr:
        log(proc.stderr.strip())

    if proc.returncode == 0 and work_out.exists():
        shutil.copy2(work_out, out_path)
        log(f"Done: {out_path}")
        try:
            summary = create_analysis_outputs(out_path)
        except Exception as exc:
            messagebox.showerror("Excel", f"AIRR TSV was created, but Excel/counts creation failed:\n{exc}")
            return
        log(f"Counts TSV: {summary['counts_tsv']}")
        log(f"Excel: {summary['xlsx']}")
        log(f"Rows included in counts: {summary['included_rows']:,}")
        log(f"Unique clonotypes: {summary['unique_clonotypes']:,}")
        messagebox.showinfo(
            "Complete",
            "Outputs:\n"
            f"{out_path}\n"
            f"{summary['counts_tsv']}\n"
            f"{summary['xlsx']}",
        )
        shutil.rmtree(work_dir, ignore_errors=True)
    else:
        log(f"IgBLAST work folder kept for troubleshooting: {work_dir}")
        messagebox.showerror("Error", "IgBLAST failed. Check the log in this window.")


def make_excel_from_existing_tsv(input_path, output_dir, log):
    if not input_path:
        messagebox.showwarning("Input", "Please select the original FASTA file first.")
        return
    in_path = Path(input_path)
    out_dir = Path(output_dir) if output_dir else RESULT_DIR
    airr_path = out_dir / (output_base_name(in_path) + ".igblast.airr.tsv")
    if not airr_path.exists():
        legacy_path = out_dir / (in_path.stem + ".igblast.airr.tsv")
        if legacy_path.exists():
            airr_path = legacy_path
        else:
            candidates = list(out_dir.glob("*.igblast.airr.tsv"))
            if len(candidates) == 1:
                airr_path = candidates[0]
            else:
                messagebox.showerror("AIRR TSV", f"Expected AIRR TSV not found:\n{airr_path}")
                return
    try:
        summary = create_analysis_outputs(airr_path)
    except Exception as exc:
        messagebox.showerror("Excel", f"Failed to create Excel/counts:\n{exc}")
        return
    log(f"Counts TSV: {summary['counts_tsv']}")
    log(f"Excel: {summary['xlsx']}")
    log(f"Rows included in counts: {summary['included_rows']:,}")
    log(f"Unique clonotypes: {summary['unique_clonotypes']:,}")
    messagebox.showinfo("Complete", f"Excel:\n{summary['xlsx']}")


def show_setup_status(log):
    missing = collect_missing_prereqs()
    if missing:
        log("Missing setup files:")
        for path in missing:
            log(path)
        messagebox.showerror("Missing files", "Missing files:\n" + "\n".join(missing))
        return
    log(f"Setup OK. Reference dir: {get_ref_dir()}")
    log(f"IgBLAST: {IGBLAST}")
    messagebox.showinfo("Setup OK", "IgBLAST and reference data were found.")


def main():
    root = tk.Tk()
    root.title("RG merged FASTA -> AIRR TSV / Excel")
    root.geometry("860x440")

    frame = tk.Frame(root, padx=10, pady=10)
    frame.pack(fill=tk.BOTH, expand=True)

    tk.Label(frame, text="Merged FASTA (PRESTO assemble-pass / extendedFrags):").grid(row=0, column=0, sticky="w")
    input_var = tk.StringVar()
    entry = tk.Entry(frame, textvariable=input_var, width=80)
    entry.grid(row=1, column=0, padx=(0, 8), sticky="we")

    output_dir_var = tk.StringVar(value=str(RESULT_DIR))

    def browse():
        path = filedialog.askopenfilename(
            title="Select FASTA",
            filetypes=[("FASTA", "*.fasta;*.fa;*.fna"), ("All", "*.*")],
        )
        if path:
            input_var.set(path)
            if output_dir_var.get().strip() == str(RESULT_DIR):
                output_dir_var.set(str(Path(path).parent))

    tk.Button(frame, text="Browse", command=browse).grid(row=1, column=1, sticky="e")

    tk.Label(frame, text="Output folder:").grid(row=2, column=0, sticky="w", pady=(8, 0))
    output_entry = tk.Entry(frame, textvariable=output_dir_var, width=80)
    output_entry.grid(row=3, column=0, padx=(0, 8), sticky="we")

    def browse_output_dir():
        path = filedialog.askdirectory(title="Select output folder")
        if path:
            output_dir_var.set(path)

    tk.Button(frame, text="Browse", command=browse_output_dir).grid(row=3, column=1, sticky="e")

    log_box = tk.Text(frame, height=10)
    log_box.grid(row=4, column=0, columnspan=2, pady=(10, 0), sticky="nsew")

    def log(msg):
        stamp = datetime.now().strftime("%H:%M:%S")
        log_box.insert(tk.END, f"[{stamp}] {msg}\n")
        log_box.see(tk.END)

    def run():
        run_igblast(input_var.get().strip(), output_dir_var.get().strip(), log)

    buttons = tk.Frame(frame)
    buttons.grid(row=5, column=0, columnspan=2, pady=10, sticky="w")
    tk.Button(buttons, text="Run IgBLAST + Excel", command=run).pack(side=tk.LEFT)
    tk.Button(
        buttons,
        text="Excel from existing TSV",
        command=lambda: make_excel_from_existing_tsv(input_var.get().strip(), output_dir_var.get().strip(), log),
    ).pack(side=tk.LEFT, padx=(8, 0))
    tk.Button(buttons, text="Check setup", command=lambda: show_setup_status(log)).pack(side=tk.LEFT, padx=(8, 0))

    frame.columnconfigure(0, weight=1)
    frame.rowconfigure(4, weight=1)

    root.mainloop()


if __name__ == "__main__":
    main()
