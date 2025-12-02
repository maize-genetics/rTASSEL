# TasselPhenotype Class Definition

Defines the `TasselPhenotype` class, which represents phenotype data in
the TASSEL 5 framework.

## Slots

- `attrData`:

  A `tbl_df` containing attribute data for the phenotype.

- `attrSummary`:

  A `list` summarizing the attributes of the phenotype data.

- `dispData`:

  A `java_pheno_tbl` object for displaying phenotype data.

- `rData`:

  A `tbl_df` containing the phenotype data in R format.

- `jRefObj`:

  A `jobjRef` representing a reference to the Java object in TASSEL 5.

- `jMemAddress`:

  A `character` string representing the memory address of the Java
  object.

- `jClass`:

  A `character` string representing the Java class name of the object.
