import { SomeFile } from "../utils.js";

let Fastq1 = new SomeFile("reads.fastq", "default", `
@r1
TGAATGCGAACTCCGGGACGCTCAGTAATGTGACGATAGCTGAAAACTGTACGATAAACNGTACGCTGAGGGCAGAAAAAATCGTCGGGGACATTNTAAAGGCGGCGAGCGCGGCTTTTCCG
+
+"@6<:27(F&5)9)"B:%B+A-%5A?2$HCB0B+0=D<7E/<.03#!.F77@6B==?C"7>;))%;,3-$.A06+<-1/@@?,26">=?*@'0;$:;??G+:#+(A?9+10!8!?()?7C>
@r2
NTTNTGATGCGGGCTTGTGGAGTTCAGCCGATCTGACTTATGTCATTACCTATGAAATGTGAGGACGCTATGCCTGTACCAAATCCTACAATGCCGGTGAAAGGTGCCGGGATCACCCTGTGGGTTTATAAGGGGATCGGTGACCCCTACGCGAATCCGCTTTCAGACGTTGACTGGTCGCGTCTGGCAAAAGTTAAAGACCTGACGCCCGGCGAACTGACCGCTGAGNCCTATGACGACAGCTATCTCGATGATGAAGATGCAGACTGGACTGC
+
(#!!'+!$""%+(+)'%)%!+!(&++)''"#"#&#"!'!("%'""("+&%$%*%%#$%#%#!)*'(#")(($&$'&%+&#%*)*#*%*')(%+!%%*"$%"#+)$&&+)&)*+!"*)!*!("&&"*#+"&"'(%)*("'!$*!!%$&&&$!!&&"(*"$&"#&!$%'%"#)$#+%*+)!&*)+(""#!)!%*#"*)*')&")($+*%%)!*)!('(%""+%"$##"#+(('!*(($*'!"*('"+)&%#&$+('**$$&+*&!#%)')'(+(!%+
`);


export const Lessons = [
	{
		id: "lesson-loading",
		inputs: [],
	},
	{
		id: "map-1",
		title: "bowtie to map reads to a reference genome",
		description: "that maps the reads in <code>reads.fastq</code> to the reference located at <code>/example/index/lambda_virus</code>",
		usage: "--help",
		tool: "",
		inputs: [Fastq1],
		goal: "-x /example/index/lambda_virus -U reads.fastq",
	},
	{
		id: "map-2",
		title: "test",
		description: "that ...",
		usage: "--help",
		tool: "",
		inputs: [Fastq1],
		goal: "-x /example/index/lambda_virus -U reads.fastq",
	},
];
