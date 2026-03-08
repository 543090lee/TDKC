# Discrete-Distilled-Model

writing in rust

TODO: 

add an option/parameter where if you do --consensus/no-fp then when we build, we look at accession information, and based on how well it's distributed or hitting many same targets, we include it or not.
Good idea


~~make load and writing bulk, minimize write_all : loead is now fixed, maybe faster write~~

~~i might be using 2x memory when loading Accession registry since it's making id_to_name and name_to_id. but i dont need name_to_id during query time.~~

~~paired end read feature: just concatenate and make it coverage for that one long string~~

~~it probably has to be streamlined i think, for writing output file. kinda inefficient i think~~

~~SIMD vectorization~~

~~stream FASTQ not load all at once~~

~~reasoning on why to use 30% coverage, maybe compare it to kraken2 how they look at actual minimizer sampling~~


~~bro i dont know what to do when a minimizer gets two different taxid when extracting from kraken2 output~~

maybe make read threshold (k-l+1)/(read length-k + 1)

it seems like accession information like when it's 12059:7, right now our query.rs is outputting only the hit.accessions of the first 12059 hit from that long run...

JUST SAY AMBIGUOUS! THE ONES THAT HAVE A TIE! Or actually just put them in unclassified