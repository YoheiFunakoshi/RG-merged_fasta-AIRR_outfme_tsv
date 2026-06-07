# RG merged FASTA → AIRR TSV / Excel (IgBLAST GUI)

RG社データのR1/R2を前段でtrim-mergeして作成したmerged FASTAを入力し、IgBLASTでAIRR outfmt 19のTSVを出力し、共同研究者が確認しやすいcounts Excelも作成するGUIです。

このツールは **マージ済み配列を1本のqueryとしてIgBLASTに渡す解析** です。R1/R2をマージしない `rg-paired-fastq-igblast-airr-tsv` とは別系統で、complete VDJ配列に近いmerged readを主解析に使う目的です。

## できること
- マージ済みFASTAを選択してIgBLASTを実行
- AIRR outfmt 19 TSVを出力
- AIRR TSVから `V候補セット + J候補セット + junction_aa` のcounts TSVを作成
- 同じcountsをExcelファイルとして作成
- 既存のAIRR TSVからExcel/countsを再作成

## 現時点の運用方針

- **AIRR TSVが主データです。** IgBLAST outfmt 19の結果を原則そのまま残し、次の解析へ渡すファイルとして扱います。
- **Excelは確認用です。** Excel/countsでは `productive == T` とcanonical `junction_aa` などで余分な配列を落とし、見やすいクローン集計として使います。
- **デフォルト入力は標準のPRESTO `assemble-pass.fasta` です。** `filtered_no_suspect_R2_adapter_like` などが付いたFASTAは補助・感度確認用として扱います。
- R1/R2をマージしない解析結果は、merged解析で落ちるR1-only/R2-only由来クローンを確認するための補助データです。
- 生成されたFASTA/AIRR TSV/Excel/counts TSVは大きく、検体データを含むためGitHubへ置かないでください。

## アプリの場所
- 起動用ショートカット: `AIRR_igblast_app.lnk`
- 本体スクリプト: `AIRR_igblast_app.pyw`
- 標準出力先: `result_AIRR_outfmat/`
- GUIでは出力先フォルダを変更できます。FASTAを選択すると、標準ではFASTAと同じフォルダが出力先になります。
- FASTA名が長い場合でもExcelで開けるよう、出力ファイル名はサンプル名を短縮して作成します。

## 前提（インストール済み）
- IgBLAST: `C:/Program Files/NCBI/igblast-1.21.0/bin/igblastn.exe`

## 参照データの場所（デスクトップ上）
`C:/Users/Yohei Funakoshi/Desktop/IgBlast用参照データ`

このフォルダは **アプリから参照される必須データ** です。削除しないでください。  
フォルダ名/場所を変更した場合は `AIRR_igblast_app.pyw` の `REF_DIR_FULL` を更新してください。  
非ASCIIパスが原因で参照エラーになる場合、アプリが `refdata` というASCII名のジャンクションを自動作成します。

### フォルダ構成と意図
- `db/` : IMGTの参照FASTAから作ったBLASTデータベース（索引ファイル群）
- `IMGT_raw/` : IMGTから取得した元のFASTA（V/D/J 生殖系列配列）
- `internal_data/` : IgBLAST同梱の内部注釈データ（IMGT番号やFWR/CDR境界の補助）
- `optional_file/` : J遺伝子のフレーム/ CDR3終端などの補助情報（-auxiliary_data で指定）

## 再現手順（参照データの作り方）

### 1) IgBLASTのインストール
- 公式: https://ncbi.github.io/igblast/
- インストール先: `C:/Program Files/NCBI/igblast-1.21.0/`

### 2) 参照データフォルダを作成
```
C:/Users/Yohei Funakoshi/Desktop/IgBlast用参照データ
```

### 3) internal_data / optional_file をコピー
IgBLASTインストール先からコピーします。
```
C:/Program Files/NCBI/igblast-1.21.0/internal_data
C:/Program Files/NCBI/igblast-1.21.0/optional_file
```

### 4) IMGTからヒトIgHのV/D/J FASTAを取得
IMGTのページから Human の IGHV / IGHD / IGHJ をダウンロードし、以下に保存します。
```
C:/Users/Yohei Funakoshi/Desktop/IgBlast用参照データ/IMGT_raw
```
ファイル名例:
- `IMGT_IGHV.fasta`
- `IMGT_IGHD.fasta`
- `IMGT_IGHJ.fasta`

### 5) IMGTヘッダを簡略化（imgt.fasta作成）
IMGTのヘッダをIgBLASTが扱いやすい形式にします。

**Python例（Perlが無い場合）**
```python
from pathlib import Path

ref = Path(r"C:/Users/Yohei Funakoshi/Desktop/IgBlast用参照データ")
inputs = [
    ref/"IMGT_raw"/"IMGT_IGHV.fasta",
    ref/"IMGT_raw"/"IMGT_IGHD.fasta",
    ref/"IMGT_raw"/"IMGT_IGHJ.fasta",
]

for inp in inputs:
    outp = inp.with_suffix(".imgt.fasta")
    with inp.open("r", encoding="utf-8", errors="ignore") as fin, outp.open("w", encoding="ascii", errors="ignore") as fout:
        for line in fin:
            if line.startswith(">"):
                parts = line[1:].strip().split("|")
                gene = parts[1].strip() if len(parts) >= 2 and parts[1].strip() else line[1:].strip()
                fout.write(">" + gene + "\n")
            else:
                seq = line.strip().upper()
                if seq:
                    fout.write(seq + "\n")
```

出力例:
- `IMGT_IGHV.imgt.fasta`
- `IMGT_IGHD.imgt.fasta`
- `IMGT_IGHJ.imgt.fasta`

### 6) makeblastdbでDB作成
```
"C:/Program Files/NCBI/igblast-1.21.0/bin/makeblastdb.exe" -parse_seqids -dbtype nucl -in IMGT_IGHV.imgt.fasta -out db/IMGT_IGHV.imgt
"C:/Program Files/NCBI/igblast-1.21.0/bin/makeblastdb.exe" -parse_seqids -dbtype nucl -in IMGT_IGHD.imgt.fasta -out db/IMGT_IGHD.imgt
"C:/Program Files/NCBI/igblast-1.21.0/bin/makeblastdb.exe" -parse_seqids -dbtype nucl -in IMGT_IGHJ.imgt.fasta -out db/IMGT_IGHJ.imgt
```

## アプリの使い方
1. `AIRR_igblast_app.lnk` をダブルクリック
2. マージ済みFASTA（PRESTO `assemble-pass.fasta` / `extendedFrags.fasta` など）を選択
3. 出力先フォルダを確認
4. `Run IgBLAST + Excel` を押す
5. 出力先フォルダに以下が作成される

```text
<sample>.igblast.airr.tsv
<sample>.airr_counts.tsv
<sample>.airr_counts.xlsx
```

既に `<sample>.igblast.airr.tsv` がある場合は、同じFASTAと出力先を指定して `Excel from existing TSV` を押すと、IgBLASTを再実行せずにcounts TSVとExcelだけを作り直せます。

## counts Excelの集計ルール

merged FASTA由来の各AIRR行について、以下を満たす行だけをcounts対象にします。

- `locus` が空欄または `IGH`
- `v_call` からアリル番号を外したV候補セットが作れる
- `j_call` からアリル番号を外したJ候補セットが作れる
- `junction_aa` がある
- `junction_aa` に stop `*` がない
- `junction_aa` が `C` で始まる
- `junction_aa` が `W` または `F` で終わる
- `junction_aa` の長さが5-40 amino acids

`productive == T` をcounts対象の必須条件にします。merged FASTAでは1本の配列に対するIgBLAST判定なので、前回の非マージR1/R2統合版よりもproductiveを採用条件に使いやすいためです。

Excelには以下のシートを作成します。

- `Summary`: 入力AIRR TSV、総行数、counts対象行数、unique clonotype数など
- `Counts`: `unique_v_gene_set + unique_j_gene_set + junction_aa` ごとの集計
- `Excluded`: counts対象から外れた理由別の行数

## 補足
- 参照DBは検体固有ではなく「ヒトIgH用の一般的DB」です。
- 参照DBを更新したい場合は、IMGTのFASTAを更新し、makeblastdbを再実行してください。
- 研究用FASTA、FASTQ、AIRR TSV、Excel出力はGitHubへ置かない前提です。

## 申し送り・今後の検討点

- GUI本体は `AIRR_igblast_app.pyw` です。Windowsの日本語パス/長いパスでIgBLASTやExcelが失敗しないよう、IgBLAST実行時は `%LOCALAPPDATA%\RGMergedFastaIgblastAirrTsv\work\...` のASCII作業フォルダを使います。
- 出力名は長すぎるFASTA名を短縮し、`<sample>.igblast.airr.tsv`, `<sample>.airr_counts.tsv`, `<sample>.airr_counts.xlsx` を作成します。
- Excel/countsのユニーク定義は `unique_v_gene_set + unique_j_gene_set + junction_aa` です。D callとC領域はユニーク定義には使いません。
- 旧ルールの `old_no_productive_filter` 出力は正式ルールではありません。混乱を避けるため通常解析には使わないでください。
- 追加すると有用な改善は、実行ごとの `run_log.txt` 自動保存、IgBLASTバージョン/参照DB作成日/入力FASTAのMD5記録、標準FASTAとfiltered補助FASTAの比較サマリです。
