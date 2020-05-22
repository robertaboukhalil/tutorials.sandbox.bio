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
		id: "loading",
		title: "Loading...",
		description: "",
		usage: null,
		tool: "",
		command: "Loading...",
		hint: "",
		inputs: [],
		goal: {
			name: ""
		}
	},
	{
		id: "intersect-1",
		title: "bedtools intersect to do xyz",
		description: "that generates overlapping regions between <code>a.bed</code> and <code>b.bed</code>",
		usage: "intersect --help",
		tool: "intersect",
		command: "bedtools intersect",
		hint: "",
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
		id: "merge-1",
		title: "bedtools merge to do abc",
		description: "that merges overlapping regions between the two <code>.bed</code> files",
		usage: "merge --help",
		tool: "merge",
		command: "bedtools merge",
		hint: "",
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
