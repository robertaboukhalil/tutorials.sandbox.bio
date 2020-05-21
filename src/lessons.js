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
		usage: null,
		tool: "",
		command: "Loading...",
		inputs: [],
		goal: {
			name: ""
		}
	},
	{
		title: "bedtools intersect to do xyz",
		description: "that returns the common regions between the files <code>a.bed</code> and <code>b.bed</code>",
		usage: "intersect --help",
		tool: "intersect",
		command: `bedtools intersect`,
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
	},
	{
		title: "bedtools merge to do abc",
		description: "that merges overlapping regions between the two <code>.bed</code> files",
		usage: "merge --help",
		tool: "merge",
		command: `bedtools merge`,
		inputs: [
			{
				name: "a.bed",
				contents: Bed1
			},
		],
		goal: {
			name: "Goal",
			contents: Bed1_intersect_Bed2
		}
	}
];
