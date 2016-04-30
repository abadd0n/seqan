//
// Created by abadd0n on 4/30/16.
//

#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>


using namespace seqan;

int main() {

    SeqFileIn sfi1(toCString(getAbsolutePath("apps/LAGAN/bovine_adenovirus_6.fa")));
    SeqFileIn sfi2(toCString(getAbsolutePath("apps/LAGAN/bovine_adenovirus_D.fa")));
    CharString meta;
    Dna5String d5s1;
    Dna5String d5s2;
    readRecord(meta, d5s1, sfi1);
    readRecord(meta, d5s2, sfi2);

    std::cout << d5s1 << std::endl;
    std::cout << d5s2 << std::endl;

    return 0;
}