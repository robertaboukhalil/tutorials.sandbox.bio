export let Tools = [
	{ id: "bedtools", description: "Wrangle .bed genomic interval files" },
	{ id: "bowtie", description: "Aligns reads to the genome" },
];


export class SomeFile
{
	constructor(name, type, contents)
	{
		this.name = name || "";
		this.type = type || "default";
		this.contents = (contents || "").trim();
    }
}
