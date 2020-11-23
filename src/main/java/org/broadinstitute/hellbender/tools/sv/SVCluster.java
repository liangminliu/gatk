package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Clusters SVs with similar breakpoints based on coordinates and supporting evidence.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Unclustered structural variants from
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCluster
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Clusters structural variants",
        oneLineSummary = "Clusters structural variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class SVCluster extends VariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String DEPTH_ONLY_INCLUDE_INTERVAL_OVERLAP_LONG_NAME = "depth-include-overlap";
    public static final String MIN_DEPTH_SIZE_LONG_NAME = "min-depth-size";

    @Argument(
            doc = "Split reads file",
            fullName = SPLIT_READ_LONG_NAME
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs file",
            fullName = DISCORDANT_PAIRS_LONG_NAME
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Sample coverage tsv",
            fullName = SAMPLE_COVERAGE_LONG_NAME
    )
    private GATKPath sampleCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 50;

    @Argument(
            doc = "Depth-only call min included intervals overlap",
            fullName = DEPTH_ONLY_INCLUDE_INTERVAL_OVERLAP_LONG_NAME,
            minValue = 0,
            maxValue = 1,
            optional = true
    )
    private double minDepthOnlyIncludeOverlap = 0.5;

    @Argument(
            doc = "Minimum depth-only call size to emit",
            fullName = MIN_DEPTH_SIZE_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minDepthOnlySize = 5000;

    private SAMSequenceDictionary dictionary;

    private final Map<String,IntervalTree<Object>> includedIntervalTreeMap = new HashMap<>();
    private VariantContextWriter writer;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;

    private SVDepthOnlyCallDefragmenter defragmenter;
    private List<SVCallRecord> nonDepthRawCallsBuffer;
    private SVClusterEngine clusterEngine;
    private BreakpointRefiner breakpointRefiner;
    private SVEvidenceCollector evidenceCollector;
    private Map<String,Double> sampleCoverageMap;
    private List<String> samplesList;
    private String currentContig;

    public static String STRANDS_ATTRIBUTE = "STRANDS";
    public static String ALGORITHMS_ATTRIBUTE = "ALGORITHMS";
    public static String START_SPLIT_READ_COUNT_ATTRIBUTE = "SR1";
    public static String END_SPLIT_READ_COUNT_ATTRIBUTE = "SR2";
    public static String DISCORDANT_PAIR_COUNT_ATTRIBUTE = "PE";
    public static String RAW_CALL_ATTRIBUTE = "RC";

    public static int RAW_CALL_ATTRIBUTE_TRUE = 1;
    public static int RAW_CALL_ATTRIBUTE_FALSE = 0;

    public static String DEPTH_ALGORITHM = "depth";

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public boolean ignoresUserIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        loadSampleCoverage();
        initializeSplitReadEvidenceDataSource();
        initializeDiscordantPairDataSource();
        loadIntervalTree();

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary);
        nonDepthRawCallsBuffer = new ArrayList<>();
        clusterEngine = new SVClusterEngineNoCNV(dictionary, false, SVClusterEngine.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END);
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap);
        evidenceCollector = new SVEvidenceCollector(splitReadSource, discordantPairSource, dictionary, null);

        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        if (!defragmenter.isEmpty() || !nonDepthRawCallsBuffer.isEmpty()) {
            processClusters();
        }
        writer.close();
        return null;
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile.toString(),
                "splitReadsFile",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void initializeDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile.toString(),
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadSampleCoverage() {
        final String fileString = sampleCoverageFile.toString();
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
            samplesList = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t")[0])
                    .collect(Collectors.toList());
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(fileString, e);
        }
    }

    private void loadIntervalTree() {
        if (getRequestedIntervals() == null) {
            for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
                includedIntervalTreeMap.put(sequence.getSequenceName(), new IntervalTree<>());
                includedIntervalTreeMap.get(sequence.getSequenceName()).put(1, sequence.getSequenceLength(), null);
            }
        } else {
            for (final SimpleInterval interval : getRequestedIntervals()) {
                includedIntervalTreeMap.putIfAbsent(interval.getContig(), new IntervalTree<>());
                includedIntervalTreeMap.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), null);
            }
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord call = SVCallRecord.create(variant);
        if (!isValidSize(call, minEventSize) || !intervalIsIncluded(call, includedIntervalTreeMap, minDepthOnlyIncludeOverlap)) {
            return;
        }

        if (!call.getContig().equals(currentContig)) {
            if (currentContig != null) {
                processClusters();
            }
            currentContig = call.getContig();
        }
        if (SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)) {
            defragmenter.add(new SVCallRecordWithEvidence(call));
        } else {
            nonDepthRawCallsBuffer.add(call);
        }
    }

    private void processClusters() {
        logger.info("Clustering contig " + currentContig);
        final Stream<SVCallRecordWithEvidence> defragmentedStream = defragmenter.getOutput().stream();
        final Stream<SVCallRecordWithEvidence> nonDepthStream = nonDepthRawCallsBuffer.stream()
                .map(SVCallRecordWithEvidence::new)
                .flatMap(this::convertInversionsToBreakends);
        //Combine and sort depth and non-depth calls because they must be added in dictionary order
        Stream.concat(defragmentedStream, nonDepthStream)
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .forEachOrdered(clusterEngine::add);
        nonDepthRawCallsBuffer.clear();
        final List<SVCallRecordWithEvidence> clusteredCalls = clusterEngine.getOutput();
        final List<SVCallRecordWithEvidence> callsWithEvidence = evidenceCollector.collectEvidence(clusteredCalls);
        final List<SVCallRecordWithEvidence> refinedCalls = callsWithEvidence.stream()
                .filter(call -> !SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call) || call.getLength() >= minDepthOnlySize)
                .map(breakpointRefiner::refineCall)
                .collect(Collectors.toList());
        final List<SVCallRecordWithEvidence> finalCalls = clusterEngine.deduplicateItems(refinedCalls);
        write(finalCalls);
        logger.info("Contig " + currentContig + " successfully clustered");
    }

    private void write(final List<SVCallRecordWithEvidence> calls) {
        calls.stream()
                .sorted(Comparator.comparing(c -> c.getPosition1AsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
    }

    private void writeVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samplesList);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        header.addMetaDataLine(new VCFFormatHeaderLine(START_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(END_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
        writer.writeHeader(header);
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        return SVCallRecordUtils.getVariantBuilder(call, samplesList).make();
    }

    private Stream<SVCallRecordWithEvidence> convertInversionsToBreakends(final SVCallRecordWithEvidence call) {
        if (!call.getType().equals(StructuralVariantType.INV)) {
            return Stream.of(call);
        }
        final SVCallRecordWithEvidence positiveBreakend = new SVCallRecordWithEvidence(call.getId(), call.getContig(),
                call.getStart(), call.getEnd(), true, true, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes(), call.getStartSplitReadSites(), call.getEndSplitReadSites(),
                call.getDiscordantPairs());
        final SVCallRecordWithEvidence negativeBreakend = new SVCallRecordWithEvidence(call.getId(), call.getContig(),
                call.getStart(), call.getEnd(), false, false, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes(), call.getStartSplitReadSites(), call.getEndSplitReadSites(),
                call.getDiscordantPairs());
        return Stream.of(positiveBreakend, negativeBreakend);
    }

    public static boolean isValidSize(final SVCallRecord call, final int minEventSize) {
        return call.getType().equals(StructuralVariantType.BND) || call.getLength() >= minEventSize;
    }

    public static <T> boolean intervalIsIncluded(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap,
                                                 final double minDepthOnlyIncludeOverlap) {
        if (SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)) {
            return intervalIsIncludedDepthOnly(call, includedIntervalTreeMap, minDepthOnlyIncludeOverlap);
        }
        return intervalIsIncludedNonDepthOnly(call, includedIntervalTreeMap);
    }

    private static <T> boolean intervalIsIncludedNonDepthOnly(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap) {
        final IntervalTree<T> startTree = includedIntervalTreeMap.get(call.getContig());
        if (startTree == null) {
            return false;
        }
        final IntervalTree<T> endTree = includedIntervalTreeMap.get(call.getContig2());
        if (endTree == null) {
            return false;
        }
        return startTree.overlappers(call.getStart(), call.getStart() + 1).hasNext()
                && endTree.overlappers(call.getEnd(), call.getEnd() + 1).hasNext();
    }

    private static <T> boolean intervalIsIncludedDepthOnly(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap,
                                                           final double minDepthOnlyIncludeOverlap) {
        final IntervalTree<T> tree = includedIntervalTreeMap.get(call.getContig());
        if (tree == null) {
            return false;
        }
        final double overlapFraction = totalOverlap(call.getStart(), call.getEnd(), tree) / (double) call.getLength();
        return overlapFraction >= minDepthOnlyIncludeOverlap;
    }

    private static <T> long totalOverlap(final int start, final int end, final IntervalTree<T> tree) {
        final Iterator<IntervalTree.Node<T>> iter = tree.overlappers(start, end);
        long overlap = 0;
        while (iter.hasNext()) {
            final IntervalTree.Node<T> node = iter.next();
            overlap += intersectionLength(start, end, node.getStart(), node.getEnd());
        }
        return overlap;
    }

    private static long intersectionLength(final int start1, final int end1, final int start2, final int end2) {
        return Math.max(0, Math.min(end1, end2) - Math.max(start1, start2) + 1);
    }
}
