//
// Created by Sebastian Proft, Anton Komissarov on 4/30/16.
//

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;

template <typename TSeedChain>
void updateSeedPositions(TSeedChain & seedChain, unsigned globalPosH, unsigned globalPosV)
{
    typedef typename Iterator<TSeedChain>::Type            TSeedChainIter;
    for (TSeedChainIter it = begin(seedChain); it != end(seedChain); ++it)
    {
        // pos
        setBeginPositionH(*it, beginPositionH(*it) + globalPosH);
        setEndPositionH(*it, endPositionH(*it) + globalPosH);
        setBeginPositionV(*it, beginPositionV(*it) + globalPosV);
        setEndPositionV(*it, endPositionV(*it) + globalPosV);

        //diag
        setUpperDiagonal(*it, upperDiagonal(*it) + globalPosH-globalPosV);
        setLowerDiagonal(*it, lowerDiagonal(*it) + globalPosH-globalPosV);
    }
}

template <typename TSeedChain, typename TSeqString, typename TSeedParams, typename TSeed>
inline bool createSeedChain( TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, TSeedParams const & seedParams, TSeed const & /*tag*/ )
{
    if ( length( seqH ) < seedParams[0] || length( seqV ) < seedParams[0] )
    {
        return false;
    }

    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;
    typedef Infix<Dna5String>::Type TInfix;

    TIndex index( seqH );
    resize( indexShape( index ), seedParams[0] );
    TInfix kmer;

    //create seedSet
    typedef SeedSet<TSeed>      TSeedSet;

    TSeedSet seedSet;
    Score<int, Simple> scoringScheme( 2, -1, -2 ); // match, mismatch, gap
    for ( unsigned qPos = 0; qPos < length( seqV ) - seedParams[0] + 1; ++qPos )
    {
        kmer = infix( seqV, qPos, qPos + seedParams[0] );
        hash(indexShape(index), begin(kmer));
        for ( auto occurence : getOccurrences( index, indexShape( index ) ) )
        {
            // try to add seed to seed set using CHAOS chaining method
            if ( !addSeed( seedSet, TSeed( occurence, qPos, seedParams[0] ), seedParams[1] /*max diag dist*/, seedParams[2] /*band width*/, scoringScheme, seqH, seqV, Chaos() ) )
                addSeed( seedSet, TSeed( occurence, qPos, seedParams[0] ), Single() ); // just add seed if CHAOS fails
        }
    }

    if ( empty (seedSet) )
        return false;

    chainSeedsGlobally( seedChain, seedSet, SparseChaining() );

    return true;
}

struct LaganOptions
{
    CharString path2file1;
    CharString path2file2;
    String<Tuple<unsigned, 3> > seedParams;
};

ArgumentParser::ParseResult
parseCommandLine(LaganOptions & options, int argc, char ** argv)
{
    // Argument parser
    ArgumentParser parser("LAGAN");

    // Set short description, version, and date.
    setShortDescription(parser, "LAGAN Algorithm");
    setVersion(parser, "1.0");
    setDate(parser, "May 2016");

    // Define usage line and long description.
    addUsageLine(parser,
                 "\"\\fIFILE_ONE\\fP\" \"\\fIFILE_TWO\\fP\" [\\fIOPTIONS\\fP]");
    addDescription(parser,
                   "This program uses LAGAN algorithm to aligns two genomes");

    // We require one argument (fa or fq file to read sequences from).
    addArgument( parser, ArgParseArgument( ArgParseArgument::STRING, "PATH" ) );

    // User can also provide second fa or fq file in case he has genomes in separate files
    // second file is arg
    addArgument( parser, ArgParseArgument( ArgParseArgument::STRING, "PATH" ) );

    addOption( parser, ArgParseOption(
            "s", "seedparams", "qGram size and its respective CHAOS parameters",
            ArgParseArgument::STRING, "PARAMS", true ) );

    addTextSection(parser, "Examples");

    // Add example
    addListItem(parser, "\\fBLAGAN\\fP \\fIpath2firstFastaFile\\fP \\fIpath2firstFastaFile\\fP",
                "Align genome from \\fIfirstFastaFile\\fP against the one from \\fIsecondFastaFile\\fP, doing one iteration with default parameters. "
                        "In this iteration the program will build seed from qGrams with size of 31 and try to chain them with diagonal distance of 5 and band width of 10.");

    // Add example
    addListItem(parser, "\\fBLAGAN\\fP \\fIpath2firstFastaFile\\fP \\fIpath2firstFastaFile\\fP \\fB-s\\fP \\fI31,5,10\\fP \\fB-s\\fP \\fI20,3,8\\fP",
                "Align genome from \\fIfirstFastaFile\\fP against the one from \\fIsecondFastaFile\\fP, doing two iterations. "
                        "In the first iteration the program will build seed from qGrams with size of 31 and try to chain them with diagonal distance of 5 and band width of 10. "
                        "In the second iteration the program will build seed from qGrams with size of 20 and try to chain them with diagonal distance  of 3 and band width of 8.");

    // Parse command line.
    ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract seed parameters
    if (isSet(parser, "seedparams"))
    {
        for (CharString paramsString : getOptionValues(parser, "seedparams"))
        {
            StringSet<CharString> params;
            strSplit( params, paramsString, EqualsChar<','>() );
            unsigned qsize,
                     chaosDiagDist,
                     chaosBandWidth;

            // parse qGram size
            if ( !lexicalCast( qsize, params[0] ) )
            {
                std::cerr << "ERROR: Invalid qGramSize " << params[0] << std::endl;
                return ArgumentParser::PARSE_ERROR;
            }
            // parse diagonal distance for CHAOS chaining
            if ( !lexicalCast(chaosDiagDist, params[1]) )
            {
                std::cerr << "ERROR: Invalid maximal diagonal distance for CHAOS chaining " << params[1] << std::endl;
                return ArgumentParser::PARSE_ERROR;
            }
            // parse band width for CHAOS chaining
            if ( !lexicalCast(chaosBandWidth, params[2]) )
            {
                std::cerr << "ERROR: Invalid band width for CHAOS chaining " << params[2] << std::endl;
                return ArgumentParser::PARSE_ERROR;
            }
            appendValue(options.seedParams, Tuple<unsigned, 3>{qsize, chaosDiagDist, chaosBandWidth} );
        }
    } else
    {
        // default values
        appendValue(options.seedParams, Tuple<unsigned, 3>{31, 5, 10} );
    }

    getArgumentValue(options.path2file1, parser, 0);
    getArgumentValue(options.path2file2, parser, 1);

    return ArgumentParser::PARSE_OK;
}

int main(int argc, char ** argv) {

    // Parse the command line.
    LaganOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // read sequences from fasta files
    // create file objects for the fasta files
    CharString filePath1 = options.path2file1;
    CharString filePath2 = options.path2file2;
    if (options.path2file1[0] != '/')
    {
        filePath1 = getAbsolutePath(toCString(filePath1));
    }
    if (options.path2file2[0] != '/')
    {
        filePath2 = getAbsolutePath(toCString(filePath2));
    }
    SeqFileIn sfi1;
    if (!open(sfi1, toCString( filePath1 ) ) )
    {
        std::cerr << "ERROR: Could not open the first file" << std::endl;
        return 1;
    }
    SeqFileIn sfi2;
    if (!open(sfi2, toCString( filePath2 ) ) )
    {
        std::cerr << "ERROR: Could not open the second file" << std::endl;
        return 1;
    }

    // read the first sequences from each file object
    StringSet<CharString> metas;
    StringSet<Dna5String> seqHs;
    StringSet<Dna5String> seqVs;
    readRecords( metas, seqHs, sfi1 );
    readRecords( metas, seqVs, sfi2 );

    Dna5String seqH = concat( seqHs, "", true );
    Dna5String seqV = concat( seqVs, "", true );

    // in this string of seeds we will store the seed chains we will find
    SeedSet<Seed<Simple>> globalSeedSet;

    // define iterator for string of seeds
    typedef Iterator<SeedSet<Seed<Simple> > >::Type TSeedSetIter;

    std::cout << "Length seqH:\t" << length(seqH) << "\tseqV:\t" << length(seqV) << std::endl;

    // we seek for seeds several times. Every next iteration happens in windows where we haven't found any seeds yet.
    // we store those seeds found in repetitive iterations in localSeedChain
    String<Seed<Simple>> localSeedChain;
    Score<int, Simple> scoringScheme( 2,-1,-2 ); // match, mismatch, gap

    for ( unsigned curParam = 0; curParam < length(options.seedParams); ++curParam )
    {
        unsigned globalGapCounter = 0;  // we count how many gaps appear in the original globalSeedSet
        unsigned pos = 0;               // current position in globalSeedSet
        unsigned infixVBegin;           // begin of the current gap in seqV
        unsigned infixHBegin;           // begin of the current gap in seqH
        unsigned infixVEnd = 0;         // end of the current gap in seqV
        unsigned infixHEnd = 0;         // end of the current gap in seqH
        unsigned nextInfixVBegin = 0;   // begin of the next gap in seqV
        unsigned nextInfixHBegin = 0;   // begin of the next gap in seqH

        Dna5String seqHInfix;
        Dna5String seqVInfix;

        unsigned const seedSetLen = length(globalSeedSet);

        for ( TSeedSetIter globalIter = begin( globalSeedSet, Standard() ); globalGapCounter < seedSetLen + 1; ++globalGapCounter, ++pos )
        {
            clear(localSeedChain);

            infixHBegin = nextInfixHBegin;
            infixVBegin = nextInfixVBegin;

            globalIter = begin( globalSeedSet, Standard() ) + pos;

            // create infixes to search for seeds
            if ( globalGapCounter == seedSetLen )
            {
                infixHEnd = length(seqH);
                infixVEnd = length(seqV);

                if ( globalGapCounter == 0 )
                {
                    seqHInfix = seqH;
                    seqVInfix = seqV;
                } else
                {
                    seqHInfix = infix(seqH, infixHBegin, infixHEnd);
                    seqVInfix = infix(seqV, infixVBegin, infixVEnd);
                }

            } else
            {
                infixHEnd = beginPositionH(*globalIter);
                infixVEnd = beginPositionV(*globalIter);
                nextInfixHBegin = endPositionH(*globalIter);
                nextInfixVBegin = endPositionV(*globalIter);
                seqHInfix = infix(seqH, infixHBegin, infixHEnd);
                seqVInfix = infix(seqV, infixVBegin, infixVEnd);
            }

            // look for seeds in infixes ..
            if (createSeedChain(localSeedChain,
                                seqHInfix,
                                seqVInfix,
                                options.seedParams[curParam],
                                Seed<Simple>() ) )
            { // .. and and seed to global seed set if any found

                // use CHAOS parameters for previous step when adding seeds to global seed set
                unsigned prevParam = curParam - 1;
                if ( curParam == 0 )
                {
                    prevParam = 0;
                }

                // update seed positions according to the position of prior global seed
                updateSeedPositions(localSeedChain, infixHBegin, infixVBegin);

                // add seeds to global seed set
                for ( Iterator<String<Seed<Simple> > >::Type localIter = begin( localSeedChain, Standard() ); localIter != end( localSeedChain, Standard() ); ++localIter )
                {
                    if ( !addSeed( globalSeedSet, *localIter, options.seedParams[prevParam][1] /*max diag dist*/, options.seedParams[prevParam][2] /*band width*/, scoringScheme, seqH, seqV, Chaos() ) )
                    {
                        addSeed( globalSeedSet, *localIter, Single() ); // just add seed if CHAOS fails
                        ++pos;
                    }
                }
            }
        }
    }

    std::cout << "Searching for seeds done" << std::endl;

    // chain seeds before doing alignment
    String<Seed<Simple> > seedChain;
    chainSeedsGlobally( seedChain, globalSeedSet, SparseChaining() );

    // define scoring schemens for alignment
    Score<int, Simple> scoringSchemeAnchor( 0, -1, -1 );
    Score<int, Simple> scoringSchemeGap( 2, -1, -1, -2 );

    // initiate alignment variable
    Align<Dna5String, ArrayGaps> alignment;
    resize( rows( alignment ), 2);
    assignSource( row( alignment, 0 ), seqH );
    assignSource( row( alignment, 1 ), seqV );
    AlignConfig<false, false, false, false> alignConfig;  // ordinary global alignment

    // align
    int result = bandedChainAlignment( alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2 );

    // print out
    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;

    return 0;
}