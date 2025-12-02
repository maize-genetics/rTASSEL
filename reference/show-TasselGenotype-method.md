# Display TasselGenotype Object

This method is used to display a summary of a \`TasselGenotype\` object.
It prints genotype display information, including the number of taxa,
number of sites, and memory address of the Java object.

## Usage

``` r
# S4 method for class 'TasselGenotype'
show(object)
```

## Arguments

- object:

  An object of class \`TasselGenotype\`.

## Details

The method utilizes the \`printGtDisp\` function to format and display
the genotype data. It extracts the necessary information from the
\`TasselGenotype\` object, including the display data (\`dispData\`),
the number of taxa and sites from the Java reference object
(\`jRefObj\`), and the memory address (\`jMemAddress\`).
