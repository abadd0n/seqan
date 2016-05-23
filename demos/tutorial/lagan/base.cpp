//![include]
#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;
//![include]

//![updateSeedPosition]
//![updateSeedPositionHead]
template <typename TSeedChain>
void updateSeedPositions(TSeedChain & seedChain, unsigned globalPosH, unsigned globalPosV)
//![updateSeedPositionHead]
{
    typedef typename Iterator<TSeedChain>::Type            TSeedChainIter;
    for (TSeedChainIter it = begin(seedChain); it != end(seedChain); ++it) {
        setBeginPositionH(*it, beginPositionH(*it) + globalPosH);
        setEndPositionH(*it, endPositionH(*it) + globalPosH);
        setBeginPositionV(*it, beginPositionV(*it) + globalPosV);
        setEndPositionV(*it, endPositionV(*it) + globalPosV);

        //diag
        setUpperDiagonal(*it, upperDiagonal(*it) + globalPosH-globalPosV);
        setLowerDiagonal(*it, lowerDiagonal(*it) + globalPosH-globalPosV);
    }
}
//![updateSeedPosition]
//![createSeedChain]
//![createSeedChainHead]
template <typename TSeedChain, typename TSeqString, typename TSeedParams, typename TSeed>
inline bool createSeedChain( TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, TSeedParams const & seedParams, TSeed const & /*tag*/ )
//![createSeedChainHead]
{
    const unsigned qgramsize = seedParams[0];
    //parameters for banded chain alignment
    const unsigned maxdiagdist = seedParams[1];
    const unsigned bandwidth = seedParams[2];
    
    if ( length( seqH ) < qgramsize || length( seqV ) < qgramsize ) {
        return false;
    }

    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;

    TIndex index( seqH );
    resize( indexShape( index ), qgramsize );

    typedef Infix<Dna5String>::Type TInfix;
    TInfix kmer;
    //![seedSet]
    //create seedSet
    typedef SeedSet<TSeed>      TSeedSet;
    TSeedSet seedSet;
    //![seedSet]
    //![scoringScheme]
    Score<int, Simple> scoringScheme( 2, -1, -2 ); // match, mismatch, gap
    //![scoringScheme]
    for ( unsigned qPos = 0; qPos < length( seqV ) - qgramsize + 1; ++qPos )
    {
        kmer = infix( seqV, qPos, qPos + qgramsize );
        hash(indexShape(index), begin(kmer));
        for ( auto occurence : getOccurrences( index, indexShape( index ) ) ) {
            // add seed to seed set using CHAOS chaining method
            if ( !addSeed( seedSet, TSeed( occurence, qPos, qgramsize ), maxdiagdist /*max diag dist*/, bandwidth /*band width*/, scoringScheme, seqH, seqV, Chaos() ) )
                addSeed( seedSet, TSeed( occurence, qPos, qgramsize ), Single() ); // just add seed if CHAOS fails
        }
    }

    if ( empty (seedSet) )
        return false;

    //chain seeds
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());

    return true;
}
//![createSeedChain]
//![parser]
//![sequences]
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
                 "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
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

    // Parse command line.
    ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "seedparams"))
    {
        for (CharString paramsString : getOptionValues(parser, "seedparams"))
        {
            StringSet<CharString> params;
            strSplit(params, paramsString, EqualsChar<','>());
            unsigned qsize,
                    chaosDiagDist,
                    chaosBandWidth;

            if ( !lexicalCast( qsize, params[0] ) )
            {
                std::cerr << "ERROR: Invalid qGramSize " << params[0] << std::endl;
                return ArgumentParser::PARSE_ERROR;
            }
            if ( !lexicalCast(chaosDiagDist, params[1]) )
            {
                std::cerr << "ERROR: Invalid maximal diagonal distance for CHAOS chaining " << params[1] << std::endl;
                return ArgumentParser::PARSE_ERROR;
            }
            if ( !lexicalCast(chaosBandWidth, params[2]) )
            {
                std::cerr << "ERROR: Invalid band width for CHAOS chaining " << params[2] << std::endl;
                return ArgumentParser::PARSE_ERROR;
            }
            appendValue(options.seedParams, Tuple<unsigned, 3>{qsize, chaosDiagDist, chaosBandWidth} );
        }
    } else {
        // default values
        appendValue(options.seedParams, Tuple<unsigned, 3>{31, 5, 10} );
    }

    getArgumentValue(options.path2file1, parser, 0);
    getArgumentValue(options.path2file2, parser, 1);

    return ArgumentParser::PARSE_OK;
}
//![parser]
//![main]
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
    if (options.path2file1[0] != '/') {
        filePath1 = getAbsolutePath(toCString(filePath1));
    }
    if (options.path2file2[0] != '/') {
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
    readRecords(metas, seqHs, sfi1);
    readRecords(metas, seqVs, sfi2);

    Dna5String seqH = concat(seqHs,"",true);
    Dna5String seqV = concat(seqVs,"",true);
    //![sequences]

    //![solution]
    // in this string of seeds we will store the seed chains we will find
    SeedSet<Seed<Simple>> globalSeedSet;

    // define iterator for the seed set
    typedef Iterator<SeedSet<Seed<Simple> > >::Type TSeedSetIter;

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

        unsigned const seedSetLen = length(globalSeedSet);

        for ( TSeedSetIter globalIter = begin( globalSeedSet, Standard() ); globalGapCounter < seedSetLen + 1; ++globalGapCounter, ++pos )
        {
            clear(localSeedChain);

            infixHBegin = nextInfixHBegin;
            infixVBegin = nextInfixVBegin;

            globalIter = begin( globalSeedSet, Standard() ) + pos;
            if (globalGapCounter == seedSetLen) {
                infixHEnd = length(seqH);
                infixVEnd = length(seqV);
            } else {
                infixHEnd = beginPositionH(*globalIter);
                infixVEnd = beginPositionV(*globalIter);
                nextInfixHBegin = endPositionH(*globalIter);
                nextInfixVBegin = endPositionV(*globalIter);
            }

            Dna5String seqVInfix = infix(seqV, infixVBegin, infixVEnd);
            Dna5String seqHInfix = infix(seqH, infixHBegin, infixHEnd);

            if (createSeedChain(localSeedChain,
                                seqHInfix,
                                seqVInfix,
                                options.seedParams[curParam],
                                Seed<Simple>() ) )
            {
                unsigned prevParam = curParam - 1;
                if ( curParam == 0 )
                {
                    prevParam = 0;
                }
                updateSeedPositions(localSeedChain, infixHBegin, infixVBegin);
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

    //create the global chain
    String<Seed<Simple> > seedChain;
    chainSeedsGlobally( seedChain, globalSeedSet, SparseChaining() );
    //![solution]

    //![align]
    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seqH);
    assignSource(row(alignment, 1), seqV);
    AlignConfig<false, false, false, false> alignConfig;  // ordinary global alignment

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;
    //![align]
    return 0;
}
//![main]