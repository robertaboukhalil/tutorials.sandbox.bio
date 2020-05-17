let Bed1 = `
chr1	10	20	a1	1	+
chr1	50	70	a2	2	+
chr1	80	90	a3	2	+
`.trim();

let Bed2 = `
chr1	30	40	b1	2	+
chr1	55	65	b2	2	+
chr1	85	120	b3	4	+
`.trim();

let Bed1_intersect_Bed2 = `
chr1	55	65	a2	2	+
chr1	85	90	a3	2	+
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
			contents: Bed1_intersect_Bed2
		}
	}
];
