"""Coverage contracts: custom columns, native adapters, BAM identity and resume."""
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

from test_multi_coverage import MODULE as cli
import test_mock_end_to_end as mock
import bam_inputs
from depth_schema import write_metabat_depth


class BamInputTests(unittest.TestCase):
    setUp = mock.MockEndToEndTests.setUp
    tearDown = mock.MockEndToEndTests.tearDown
    base_run = mock.MockEndToEndTests.base_run
    run_cli = mock.MockEndToEndTests.run_cli

    def bam(self):
        path=self.root/"sample.bam";path.write_bytes(b"synthetic")
        Path(str(path)+".bai").write_bytes(b"index")
        return path

    def test_custom_columns_preserve_pairing_and_infer_lengths(self):
        self.depth.write_text("junk\tidentifier\tBmean\tAvar\tAmean\tBvar\nx\tc1\t2\t7\t3\t11\n")
        text,meta=cli.normalized_depth_inputs([cli.DepthInput(self.depth)],{"c1":100},"s",
             contig_column="identifier",coverage_columns=["Amean","Bmean"],variance_columns=["Avar","Bvar"])
        self.assertEqual(text.splitlines()[1],"c1\t100\t5\t3\t7\t2\t11")
        self.assertEqual(meta["variance_available"],[True,True])
        for options in [dict(contig_column="bad",coverage_columns=["Amean"]),
                        dict(contig_column="identifier",coverage_columns=["Amean","Amean"]),
                        dict(contig_column="identifier",coverage_columns=["Amean","Bmean"],variance_columns=["Avar"])]:
            with self.assertRaises(ValueError):
                cli.normalized_depth_inputs([cli.DepthInput(self.depth)],{"c1":100},"s",**options)

    def test_native_mean_and_observed_variance_formats(self):
        mean=self.root/"mean.tsv";observed=self.root/"observed.tsv"
        self.depth.write_text("contigName\tcontigLen\ttotalAvgDepth\tA\tA-var\nc1\t100\t3\t3\tNA\n")
        self.assertEqual(write_metabat_depth(self.depth,mean),"means")
        self.assertEqual(mean.read_text(),"contigName\tA\nc1\t3\n")
        self.depth.write_text(self.depth.read_text().replace("NA","0"))
        self.assertEqual(write_metabat_depth(self.depth,observed),"jgi")
        self.assertEqual(observed.read_bytes(),self.depth.read_bytes())
        self.depth.write_text(self.depth.read_text()+"c2\t100\t2\t2\tNA\n")
        with self.assertRaisesRegex(ValueError,"cannot mix"):
            write_metabat_depth(self.depth,mean)

    def test_bam_headers_lengths_md5_and_indexes(self):
        bam=self.bam();seq="ACGT"*25
        good="@HD\tSO:coordinate\n@SQ\tSN:c1\tLN:100\tM5:"+hashlib.md5(seq.encode()).hexdigest()+"\n"
        def check(header,lengths=None):
            def output(command):
                if "view" in command:return header
                if "idxstats" in command:return "c1\t100\t1\t0\n*\t0\t0\t0\n"
                return ""
            with patch.object(bam_inputs,"checked_output",side_effect=output):
                return bam_inputs.validate_bams([bam],lengths or {"c1":100},io.StringIO(">c1\n"+seq+"\n"))
        self.assertEqual(check(good)["reference_sequence_check"],"M5_verified")
        self.assertIn("provenance_required",check(good.split("\tM5:")[0]+"\n")["reference_sequence_check"])
        for header in [good.replace("coordinate","queryname"),good.replace("LN:100","LN:99"),
                       good.replace(hashlib.md5(seq.encode()).hexdigest(),"0"*32),good+"@SQ\tSN:c1\tLN:100\n"]:
            with self.assertRaises(ValueError):check(header)
        with self.assertRaisesRegex(ValueError,"ID or length mismatch"):check(good,{"absent":100})
        Path(str(bam)+".bai").unlink()
        with self.assertRaisesRegex(ValueError,"index"):check(good)

    def test_bam_dictionary_and_name_collisions(self):
        a=self.bam();b=self.root/"b.bam";b.write_bytes(b"b");Path(str(b)+".bai").write_bytes(b"index")
        def output(command):
            if "view" in command:
                return "@HD\tSO:coordinate\n@SQ\tSN:c1\tLN:100\n"+("@SQ\tSN:extra\tLN:10\n" if command[-1]==str(b.resolve()) else "")
            if "idxstats" in command:return "c1\t100\t0\t0\n"
            return ""
        with patch.object(bam_inputs,"checked_output",side_effect=output):
            with self.assertRaisesRegex(ValueError,"dictionaries differ"):
                bam_inputs.validate_bams([a,b],{"c1":100},io.StringIO(""))
        with self.assertRaisesRegex(ValueError,"unique"):
            bam_inputs.validate_bams([a,a],{},io.StringIO(""))

    def test_mutually_exclusive_inputs_and_selectors(self):
        output=self.root/"out";base=self.base_run(output)
        with self.assertRaises(SystemExit), patch("sys.stderr",io.StringIO()):
            cli.build_parser().parse_args([*base,"--bam",str(self.bam())])
        i=base.index("--depth");base[i:i+2]=["--bam",str(self.bam())]
        args=cli.build_parser().parse_args([*base,"--contig-column","ID","--coverage-columns","mean"])
        with self.assertRaisesRegex(ValueError,"only to --depth"):cli.validate_args(args)

    def test_samtools_index_fallback_is_not_validation_success(self):
        result=subprocess.CompletedProcess([],0,"c1\t100\t0\t0\n",
                         "[E::hts_idx_load3] Could not load index; reverting to slow method")
        with patch.object(bam_inputs.subprocess,"run",return_value=result):
            with self.assertRaisesRegex(ValueError,"Could not load index"):
                bam_inputs.checked_output(["samtools","idxstats","bad.bam"])

    def test_corrupt_coverage_cache_is_rejected(self):
        depth=self.root/"cached.tsv";depth.write_text("bad")
        depth.with_suffix(".sha256").write_text("0"*64)
        with self.assertRaisesRegex(ValueError,"modified"),patch.object(bam_inputs.subprocess,"run") as run:
            bam_inputs.calculate_depth([],depth,1)
        run.assert_not_called()

    def test_bam_end_to_end_and_resume_without_remapping(self):
        bam=self.bam()
        script = r"""#!/usr/bin/env python3
import sys
from pathlib import Path
name=Path(sys.argv[0]).name;args=sys.argv[1:]
if name=='samtools':
 if 'view' in args:print('@HD\tSO:coordinate\n'+''.join('@SQ\tSN:c%d\tLN:100\n'%i for i in (1,2,3)),end='')
 elif 'idxstats' in args:print(''.join('c%d\t100\t1\t0\n'%i for i in (1,2,3)),end='')
elif name=='jgi_summarize_bam_contig_depths':
 Path(args[args.index('--outputDepth')+1]).write_text('contigName\tcontigLen\ttotalAvgDepth\tsample\tsample-var\n'+''.join('c%d\t100\t3\t3\t7\n'%i for i in (1,2,3)))
elif name=='concoct_coverage_table.py':
 print('contig_id\tsample')
 for row in Path(args[0]).read_text().splitlines():print(row.split('\t')[3]+'\t9')
"""
        for name in ["samtools","jgi_summarize_bam_contig_depths","concoct_coverage_table.py"]:
            tool=self.tools/name;tool.write_text(script);tool.chmod(0o755)
        output=self.root/"bam_run";args=self.base_run(output)
        i=args.index("--depth");args[i:i+2]=["--bam",str(bam)]
        args += ["--binners","concoct,metabat2","--qc-tools","eukcc","--eukcc-db",str(self.eukcc_db)]
        first=self.run_cli(*args,concoct_success=True,expected_coverage_columns=("sample",))
        self.assertEqual(first.returncode,0,first.stderr)
        summary=json.loads((output/"summary.json").read_text())
        self.assertEqual(summary["completed_binners"],["concoct","metabat2"])
        fragment=next(output.glob("results/*/*/results/union/concoct/fragment_coverage.tsv"))
        self.assertIn("c1\t9",fragment.read_text())
        self.assertIn("\t3\t7",(output/"work/input/depth.jgi.tsv").read_text())
        second=self.run_cli(*args,concoct_success=True,expected_coverage_columns=("sample",))
        self.assertEqual(second.returncode,0,second.stderr)
        self.assertIn("already complete",second.stderr)
        Path(str(bam)+".bai").write_bytes(b"changed-index")
        changed=self.run_cli(*args,concoct_success=True)
        self.assertNotEqual(changed.returncode,0)
        self.assertIn("inputs, parameters or software differ",changed.stderr)
