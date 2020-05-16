let Bed1 = `
chr1	10	20	a1	1	+
chr1	100	200	a2	2	-
`.trim();

let Bed2 = `
chr1	15	30	b1	1	+
chr1	90	101	b2	2	-
chr1	100	110	b3	3	+
chr1	200	210	b4	4	+
`.trim();

let Bed1_intersect_Bed2 = `
chr1	15	20	a1	1	+
chr1	100	101	a2	2	-
chr1	100	110	a2	2	-
`.trim();

export const Lessons = [
	{
		title: "Loading...",
		description: "",
		command: "Loading...",
		inputs: [],
		goal: {
			name: "Loading..."
		},
		usage: null
	},
	{
		title: "Lesson 1: bedtools intersect",
		description: "Using bedtools intersect to do xyz",
		usage: "intersect --help",
		command: `bedtools intersect -a a.bed -b b.bed`,
		inputs: [
			{
				name: "a.bed",
				contents: Bed1
			},
			{
				name: "b.bed",
				contents: Bed2
			}
		],
		goal: {
			name: "Goal",
			color: "purple",
			contents: Bed1_intersect_Bed2
		}
	}
];
