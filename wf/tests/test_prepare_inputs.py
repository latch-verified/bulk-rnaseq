# from wf.models import LatchGenome, Sample, SingleEndReads, Strandedness
# from wf.prepare_inputs import prepare_inputs, rsem_prepare_reference
# from wf.types import FastaFile, GtfFile, LatchFile
#
#
#
# # zipped/unzipped single/paired
# # with lgenome
# # with custom genome + gtf
# # with custom trans + gtf
#
#
# unzipped = [
#     Sample(
#         name="test",
#         strandedness=Strandedness.auto,
#         replicates=[
#             SingleEndReads(
#                 r1=LatchFile("/Users/kenny/latch/latch-verified/genes/genes.fa"),
#             ),
#         ],
#     ),
# ]
#
# prepare_inputs(
#     samples=[],
#     run_name="",
#     latch_genome=LatchGenome.RefSeq_hg38_p14,
#     custom_output_dir=None,
#     custom_gtf=None,
#     custom_ref_genome=None,
#     custom_ref_trans=None,
#     custom_salmon_index=None,
#     run_splicing=False,
# )
