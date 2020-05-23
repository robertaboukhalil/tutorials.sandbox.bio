export class BedFile
{
	constructor(name, contents, type)
	{
		this.name = name || "";
		this.contents = (contents || "").trim();
		this.type = type || "default";
    }
}
