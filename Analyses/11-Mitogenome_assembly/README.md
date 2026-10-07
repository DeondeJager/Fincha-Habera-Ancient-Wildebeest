# Assembly of Fossil936 and modern mitogenomes
- The `Fossil936` folder contains the paleomix `yaml`file and `slurm` script for mapping reads against the final Fossil936 assembly to validate it.
- The `Liu2024_mitogenomes` folder contains the [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) scripts used to assemble the mitogenomes.
  - The subfolder `Liu2024_scripts` contains the per-sample GetOrganelle scripts, which are run as an slurm array with `GetOrganelle_Liu2024.batch.slurm`.
  - Also in this folder is the `GetOrganelle_Liu2024.processResults.sh` script, which:
    
	i. Gets stats for all the assembled mitogenomes,
	
	ii. Annotates the mitogenomes with [MITOS2](https://gitlab.com/Bernt/MITOS/-/tree/mitos2) and gets the start position of the `trnF` gene, which is where the reference mitogenomes "start"
	
	iii. "Rolls" the mitogenomes to start from this `trnF` position using the `circules.py` script from [MITObim](https://github.com/chrishah/MITObim) (included here), so they all start at the same position as each other and the reference mitogenomes,
	
	iv. Re-annotates the rolled mitogenomes with MITOS2.
