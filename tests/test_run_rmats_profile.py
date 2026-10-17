"""Conservative rMATS is isolated from the original novelty-enabled route."""
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src import run_rmats_profile as profile
from src.umbrella_modules import module_reference_id, parse_selection


class RmatsProfileTests(unittest.TestCase):
    def test_selection_is_opt_in_and_does_not_change_existing_ids(self):
        self.assertEqual(parse_selection({})[0], [])
        self.assertEqual(parse_selection({'umbrella_modules': ['rmats_conservative']})[0],
                         ['rmats_conservative'])
        for tool in ('majiq', 'qapa', 'dapars2', 'vasttools'):
            self.assertEqual(module_reference_id('ref', tool, {}), module_reference_id(
                'ref', tool, {'umbrella_modules': ['rmats_conservative']}))

    def test_native_prep_omits_novelss_and_records_processing_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def native(command, **kwargs):
                self.assertNotIn('--novelSS', command)
                self.assertEqual(command[command.index('--libType')+1], 'fr-unstranded')
                tmp = Path(command[command.index('--tmp')+1])
                (tmp/'generated.rmats').write_text('native summary')
            with patch.object(profile.subprocess, 'run', side_effect=native):
                profile.prep(root/'aligned.bam', root/'annotation.gtf', 100, root/'cache', 'run1')
            self.assertEqual((root/'cache/run.rmats').read_text(), 'native summary')
            record=json.loads((root/'cache/cache.json').read_text())
            self.assertEqual(record['tool_family'], 'rmats')
            self.assertFalse(record['novel_splice_sites'])

    def test_normalization_preserves_native_effect_p_and_fdr(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);native=root/'SE.MATS.JC.txt'
            native.write_text('ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\t'
                'upstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tPValue\tFDR\tIncLevelDifference\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\n'
                '1\t"G1"\t"gene"\tchr1\t+\t200\t250\t50\t100\t300\t350\t0.002\t0.03\t-0.4\t10,11\t20,21\t30,31\t2,3\n')
            rows=profile.normalized_rows(native,'SE',{'replicates':{'a':['a1','a2'],'b':['b1','b2']}},{})
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['gene_id'],'G1')
            self.assertEqual((rows[0]['start'],rows[0]['end']),(201,250))
            self.assertEqual(float(rows[0]['effect']),-.4)
            self.assertEqual(float(rows[0]['p_value']),.002)
            self.assertEqual(float(rows[0]['q_value']),.03)
            self.assertEqual(rows[0]['status'],'tested')
            native.write_text(native.read_text().replace('0.002\t0.03','NA\tNA'))
            self.assertEqual(profile.normalized_rows(native,'SE',{'replicates':{'a':[],'b':[]}},{})[0]['status'],'untested')

    def test_event_identity_does_not_depend_on_individual_junction_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            native=Path(directory)/'SE.MATS.JC.txt'
            header=('ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\t'
                'upstreamEE\tdownstreamES\tdownstreamEE\tPValue\tFDR\tIncLevelDifference\t'
                'IJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tupstream_to_target_count\n')
            row='1\tG1\tgene\tchr1\t+\t200\t250\t50\t100\t300\t350\t0.002\t0.03\t-0.4\t10,11\t20,21\t30,31\t2,3\t5,6,7,8\n'
            native.write_text(header+row)
            preflight={'replicates':{'a':['a1','a2'],'b':['b1','b2']}}
            before=profile.normalized_rows(native,'SE',preflight,{})[0]['feature_id']
            native.write_text(header+row.replace('5,6,7,8','50,60,70,80'))
            self.assertEqual(before,profile.normalized_rows(native,'SE',preflight,{})[0]['feature_id'])

    def test_missing_biological_replication_is_not_pseudoreplicated(self):
        preflight={'inference_supported':True,'tools':{'rmats':{'supported':True,'reasons':[]}},
                   'collapse':{'r1':'a','r2':'a','r3':'b','r4':'c'},
                   'replicates':{'a':['a'],'b':['b','c']}}
        with tempfile.TemporaryDirectory() as directory, patch.object(profile.subprocess,'run') as native:
            profile.compare(preflight,{},'unused.gtf',directory)
            native.assert_not_called()
            status=json.loads((Path(directory)/'status.json').read_text())
            self.assertEqual(status['status'],'unsupported')
            self.assertTrue((Path(directory)/'results.tsv.gz').exists())

    def comparison_fixture(self, root):
        import hashlib
        runs=['r1','r2','r3','r4']
        preflight={'inference_supported':True,'tools':{'rmats':{'supported':True,'reasons':[]}},
                   'collapse':dict(zip(runs,runs)), 'replicates':{'a':runs[:2],'b':runs[2:]},
                   'runs':{'a':runs[:2],'b':runs[2:]},'layouts':['PE'],
                   'comparison_id':'a_vs_b','project_id':'project'}
        caches={}
        for run in runs:
            path=root/run;path.mkdir();(path/'run.rmats').write_text(run)
            record=dict(profile.POLICY,run_id=run,bam='/original/'+run+'.bam',layout='PE',
                        read_length=100,sha256=hashlib.sha256(run.encode()).hexdigest())
            (path/'cache.json').write_text(json.dumps(record));caches[run]=str(path)
        return preflight,caches

    def test_post_reuses_prep_without_opening_bam_and_keeps_native_evidence(self):
        import tarfile
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);preflight,caches=self.comparison_fixture(root)
            old=root/'original-result';old.write_text('preserve me')
            def native(command,**kwargs):
                self.assertNotIn('--novelSS',command)
                work=Path(command[command.index('--od')+1])
                group=Path(command[command.index('--b1')+1]).read_text().strip()
                self.assertEqual(group,'/original/r1.bam,/original/r2.bam')
                if command[-1]=='post':
                    for kind in profile.KINDS:
                        (work/(kind+'.MATS.JC.txt')).write_text('ID\tGeneID\n')
                    (work/'fromGTF.SE.txt').write_text('native event definition\n')
                    (work/'individualCounts.SE.txt').write_text('junction breakdown\n')
            with patch.object(profile.subprocess,'run',side_effect=native) as calls:
                profile.compare(preflight,caches,'annotation.gtf',root/'new')
            self.assertEqual(calls.call_count,2)
            with tarfile.open(root/'new/native.tar.gz') as t:
                self.assertEqual(t.extractfile('individualCounts.SE.txt').read(),b'junction breakdown\n')
            self.assertEqual(old.read_text(),'preserve me')
            self.assertEqual(json.loads((root/'new/status.json').read_text())['status'],'empty')

    def test_cache_tampering_fails_before_native_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);preflight,caches=self.comparison_fixture(root)
            (Path(caches['r1'])/'run.rmats').write_text('corrupted')
            with patch.object(profile.subprocess,'run') as native:
                with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                    profile.compare(preflight,caches,'annotation.gtf',root/'new')
                native.assert_not_called()

    def provenance_fixture(self, root):
        preflight,caches=self.comparison_fixture(root)
        ids={'reference_id':'ref1','module_reference_id':'mod1','analysis_id':'an1'}
        moved={}
        for run,path in caches.items():
            record=json.loads((Path(path)/'cache.json').read_text())
            record.update(reference_id='ref1',module_reference_id='mod1',cache_id='c'+run)
            target=root/'cache'/run/('c'+run);target.mkdir(parents=True)
            (target/'run.rmats').write_text(run)
            (target/'cache.json').write_text(json.dumps(record));moved[run]=str(target)
        return preflight,moved,ids

    def test_cache_from_another_reference_profile_or_layout_is_refused(self):
        for key,value,message in (('reference_id','ref2','reference_id ref2'),
                                  ('module_reference_id','mod2','module_reference_id mod2'),
                                  ('cache_id','cother','cache_id cother'),
                                  ('layout','SE','is SE, the comparison is PE')):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                root=Path(directory);preflight,caches,ids=self.provenance_fixture(root)
                path=Path(caches['r2'])/'cache.json'
                record=json.loads(path.read_text());record[key]=value;path.write_text(json.dumps(record))
                with patch.object(profile.subprocess,'run') as native:
                    with self.assertRaisesRegex(ValueError,message):
                        profile.compare(preflight,caches,'annotation.gtf',root/'new',ids=ids)
                    native.assert_not_called()

    def test_matching_provenance_runs_and_mixed_read_lengths_use_the_maximum(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);preflight,caches,ids=self.provenance_fixture(root)
            path=Path(caches['r3'])/'cache.json'
            record=json.loads(path.read_text());record['read_length']=150;path.write_text(json.dumps(record))
            lengths=[]
            header=('ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\t'
                    'upstreamEE\tdownstreamES\tdownstreamEE\tPValue\tFDR\tIncLevelDifference\t'
                    'IJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\n')
            def native(command,**kwargs):
                lengths.append(command[command.index('--readLength')+1])
                work=Path(command[command.index('--od')+1])
                if command[-1]=='post':
                    for kind in profile.KINDS:
                        (work/(kind+'.MATS.JC.txt')).write_text('ID\tGeneID\n' if kind!='SE' else header+
                            '1\tG1\tgene\tchr1\t+\t200\t250\t50\t100\t300\t350\t0.002\t0.03\t-0.4\t1,1\t2,2\t3,3\t4,4\n')
            with patch.object(profile.subprocess,'run',side_effect=native):
                profile.compare(preflight,caches,'annotation.gtf',root/'new',ids=ids)
            self.assertEqual(lengths,['150','150'])
            status=json.loads((root/'new/status.json').read_text())
            self.assertEqual((status['status'],status['tested']),('ok',1))
            self.assertEqual(status['prep_read_lengths'],[100,150])
            self.assertIn('mixed prep read lengths',status['warnings'][0])

    def test_events_without_any_test_report_no_tests(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);preflight,caches=self.comparison_fixture(root)
            header=('ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\t'
                    'upstreamEE\tdownstreamES\tdownstreamEE\tPValue\tFDR\tIncLevelDifference\t'
                    'IJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\n')
            def native(command,**kwargs):
                work=Path(command[command.index('--od')+1])
                if command[-1]=='post':
                    for kind in profile.KINDS:
                        (work/(kind+'.MATS.JC.txt')).write_text('ID\tGeneID\n' if kind!='SE' else header+
                            '1\tG1\tgene\tchr1\t+\t200\t250\t50\t100\t300\t350\tNA\tNA\tNA\t0,0\t0,0\t0,0\t0,0\n')
            with patch.object(profile.subprocess,'run',side_effect=native):
                profile.compare(preflight,caches,'annotation.gtf',root/'new')
            status=json.loads((root/'new/status.json').read_text())
            self.assertEqual((status['status'],status['events'],status['tested']),('no_tests',1,0))

if __name__=='__main__':
    unittest.main()
