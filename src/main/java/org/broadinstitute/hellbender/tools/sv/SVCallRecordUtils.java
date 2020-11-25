package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class SVCallRecordUtils {

    /**
     * Create a variant from a call for VCF interoperability
     *
     * @param call variant to convert
     * @return
     */
    public static VariantContextBuilder getVariantBuilder(final SVCallRecord call) {
        Utils.nonNull(call);
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final Allele refAllele = Allele.REF_N;
        final int end = isIntrachromosomal(call) ? call.getPositionB() : call.getPositionA() + 1;
        final VariantContextBuilder builder = new VariantContextBuilder(call.getId(), call.getContigA(), call.getPositionA(),
                end, Lists.newArrayList(refAllele, altAllele));
        builder.id(call.getId());
        builder.attribute(VCFConstants.END_KEY, end);
        builder.attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, call.getContigB());
        builder.attribute(GATKSVVCFConstants.END2_ATTRIBUTE, call.getPositionB());
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        builder.attribute(GATKSVVCFConstants.SVTYPE, call.getType());
        builder.attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, getStrandString(call));
        builder.attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, call.getAlgorithms());
        builder.genotypes(call.getGenotypes());
        return builder;
    }

    public static VariantContextBuilder fillMissingGenotypes(final VariantContextBuilder builder, final Set<String> samples) {
        Utils.nonNull(samples);
        Utils.containsNoNull(samples, "Null sample found");
        final Set<String> missingSamples = Sets.difference(samples, builder.getGenotypes().getSampleNames());
        if (!missingSamples.isEmpty()) {
            final List<Genotype> genotypes = new ArrayList<>(builder.getGenotypes().size() + missingSamples.size());
            genotypes.addAll(builder.getGenotypes());
            for (final String sample : missingSamples) {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                genotypes.add(genotypeBuilder.make());
            }
            builder.genotypes(genotypes);
        }
        return builder;
    }

    public static VariantContextBuilder wipeGenotypesAndSetRawCallAttribute(final VariantContextBuilder builder,
                                                                            final SVCallRecord call) {
        Utils.nonNull(builder);
        Utils.nonNull(call);
        final List<Genotype> genotypes = new ArrayList<>();
        for (final Genotype genotype : builder.getGenotypes()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype.getSampleName());
            if (SVCallRecord.isCarrier(genotype) || SVCallRecord.isRawCall(genotype)) {
                genotypeBuilder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
            } else {
                genotypeBuilder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
            }
            genotypes.add(genotypeBuilder.make());
        }
        return builder.genotypes(genotypes);
    }

    public static VariantContextBuilder getVariantWithEvidenceBuilder(final SVCallRecordWithEvidence call) {
        final VariantContextBuilder builder = getVariantBuilder(call);
        final SplitReadSite startSplitReadCounts = getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getPositionA());
        final SplitReadSite endSplitReadCounts = getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getPositionB());
        final Map<String,Integer> discordantPairCounts = getDiscordantPairCountsMap(call.getDiscordantPairs());
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size());
        final boolean isDepthOnly = SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call);
        for (final Genotype genotype : genotypes) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            if (!isDepthOnly) {
                final int startCount = startSplitReadCounts == null ? 0 : startSplitReadCounts.getCount(sample);
                final int endCount = endSplitReadCounts == null ? 0 : endSplitReadCounts.getCount(sample);
                genotypeBuilder.attribute(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, startCount);
                genotypeBuilder.attribute(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, endCount);
                genotypeBuilder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, discordantPairCounts.getOrDefault(sample, 0));
            }
            if (call.getCalledSamples().contains(sample)) {
                genotypeBuilder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
            } else {
                genotypeBuilder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
            }
            newGenotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(newGenotypes);
        return builder;
    }

    private static String getStrandString(final SVCallRecord call) {
        return getStrandString(call.getStrandA()) + getStrandString(call.getStrandB());
    }

    private static String getStrandString(final boolean strand) {
        return strand ? SVCallRecord.STRAND_PLUS : SVCallRecord.STRAND_MINUS;
    }

    private static SplitReadSite getSplitReadCountsAtPosition(final List<SplitReadSite> sites, final int pos) {
        Utils.nonNull(sites);
        Utils.validateArg(pos > 0, "Non-positive position");
        if (sites.stream().map(SplitReadSite::getPosition).distinct().count() != sites.size()) {
            throw new IllegalArgumentException("Sites did not have unique positions");
        }
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .findAny()
                .orElse(null);
    }

    private static Map<String,Integer> getDiscordantPairCountsMap(final Collection<DiscordantPairEvidence> discordantPairs) {
        Utils.nonNull(discordantPairs);
        return discordantPairs.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.reducing(0, e -> 1, Integer::sum)));
    }

    public static <T extends Locatable2D> Comparator<T> getLocatable2DComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareLocatable2D((T) o1, (T) o2, dictionary);
    }

    public static <T extends SVCallRecord> Comparator<T> getSiteInfoComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareRecordsBySiteInfo((T) o1, (T) o2, dictionary);
    }

    public static int compareLocatable2D(final Locatable2D first, final Locatable2D second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        // First locus
        final Comparator<Locatable> locatableComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final int compareA = locatableComparator.compare(new SimpleInterval(first.getContigA(), first.getPositionA(), first.getPositionA()),
                new SimpleInterval(second.getContigA(), second.getPositionA(), second.getPositionA()));
        if (compareA != 0) return compareA;
        // Second locus
        final int compareB = locatableComparator.compare(new SimpleInterval(first.getContigB(), first.getPositionB(), first.getPositionB()),
                new SimpleInterval(second.getContigB(), second.getPositionB(), second.getPositionB()));
        return compareB;
    }

    public static int compareRecordsBySiteInfo(final SVCallRecord first, final SVCallRecord second, final SAMSequenceDictionary dictionary) {
        final int compareLocatables = compareLocatable2D(first, second, dictionary);
        if (compareLocatables != 0) return compareLocatables;

        //Strands
        final int compareStartStrand = Boolean.compare(first.getStrandA(), first.getStrandA());
        if (compareStartStrand != 0) return compareStartStrand;
        final int compareEndStrand = Boolean.compare(first.getStrandB(), first.getStrandB());
        if (compareEndStrand != 0) return compareEndStrand;

        // Length
        final int compareLength = Integer.compare(first.getLength(), second.getLength());
        if (compareLength != 0) return compareLength;

        // Type
        final int compareType = first.getType().compareTo(second.getType());
        return compareType;
    }

    public static boolean isValidSize(final SVCallRecord call, final int minEventSize) {
        return call.getType().equals(StructuralVariantType.BND) || call.getLength() >= minEventSize;
    }

    public static <T> boolean intervalIsIncluded(final SVCallRecord call, final Map<String, IntervalTree<T>> includedIntervalTreeMap,
                                                 final double minDepthOnlyIncludeOverlap) {
        if (SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)) {
            return intervalIsIncludedDepthOnly(call, includedIntervalTreeMap, minDepthOnlyIncludeOverlap);
        }
        return intervalIsIncludedNonDepthOnly(call, includedIntervalTreeMap);
    }

    private static <T> boolean intervalIsIncludedNonDepthOnly(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap) {
        final IntervalTree<T> startTree = includedIntervalTreeMap.get(call.getContigA());
        if (startTree == null) {
            return false;
        }
        final IntervalTree<T> endTree = includedIntervalTreeMap.get(call.getContigB());
        if (endTree == null) {
            return false;
        }
        return startTree.overlappers(call.getPositionA(), call.getPositionA() + 1).hasNext()
                && endTree.overlappers(call.getPositionB(), call.getPositionB() + 1).hasNext();
    }

    private static <T> boolean intervalIsIncludedDepthOnly(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap,
                                                           final double minDepthOnlyIncludeOverlap) {
        final IntervalTree<T> tree = includedIntervalTreeMap.get(call.getContigA());
        if (tree == null) {
            return false;
        }
        final double overlapFraction = totalOverlap(call.getPositionA(), call.getPositionB(), tree) / (double) call.getLength();
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

    public static Stream<SVCallRecord> convertInversionsToBreakends(final SVCallRecord call) {
        if (!call.getType().equals(StructuralVariantType.INV)) {
            return Stream.of(call);
        }
        Utils.validateArg(isIntrachromosomal(call), "Inversion is not intrachromosomal");
        final SVCallRecord positiveBreakend = new SVCallRecord(call.getId(), call.getContigA(),
                call.getPositionA(), call.getPositionB(), true, true, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes());
        final SVCallRecord negativeBreakend = new SVCallRecord(call.getId(), call.getContigA(),
                call.getPositionA(), call.getPositionB(), false, false, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes());
        return Stream.of(positiveBreakend, negativeBreakend);
    }

    public static <T extends Locatable2D> List<T> deduplicateLocatables(final List<T> items,
                                                                      final SAMSequenceDictionary dictionary,
                                                                      final BiPredicate<T,T> itemsAreIdenticalFunction,
                                                                      final Function<Collection<T>,T> deduplicateIdenticalItemsFunction) {
        final List<T> sortedItems = items.stream().sorted(getLocatable2DComparator(dictionary)).collect(Collectors.toList());
        final List<T> deduplicatedList = new ArrayList<>();
        int i = 0;
        while (i < sortedItems.size()) {
            final T record = sortedItems.get(i);
            int j = i + 1;
            final Collection<Integer> identicalItemIndexes = new ArrayList<>();
            while (j < sortedItems.size() && record.getPositionA() == sortedItems.get(j).getPositionB()) {
                final T other = sortedItems.get(j);
                if (itemsAreIdenticalFunction.test(record, other)) {
                    identicalItemIndexes.add(j);
                }
                j++;
            }
            if (identicalItemIndexes.isEmpty()) {
                deduplicatedList.add(record);
                i++;
            } else {
                identicalItemIndexes.add(i);
                final List<T> identicalItems = identicalItemIndexes.stream().map(sortedItems::get).collect(Collectors.toList());
                deduplicatedList.add(deduplicateIdenticalItemsFunction.apply(identicalItems));
                i = j;
            }
        }
        return deduplicatedList;
    }

    public static void validateCoordinates(final SVCallRecord call) {
        Utils.nonNull(call);
        Utils.validateArg(call.getPositionA() >= 1, "Call start non-positive");
        if (isIntrachromosomal(call)) {
            Utils.validateArg(call.getPositionA() <= call.getPositionB(), "Second position before end on same contig");
        } else {
            Utils.validateArg(call.getPositionB() >= 1, "Call second position non-positive");
        }
    }

    public static void validateCoordinatesWithDictionary(final SVCallRecord call, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary);
        validateCoordinates(call);
        final SAMSequenceRecord contigARecord = dictionary.getSequence(call.getContigA());
        Utils.validateArg(contigARecord != null, "Call first contig " + call.getContigA() + " not in dictionary");
        final SAMSequenceRecord contigBRecord = dictionary.getSequence(call.getContigB());
        Utils.validateArg(contigBRecord != null, "Call second contig " + call.getContigB() + " not in dictionary");
        Utils.validateArg(call.getPositionA() <= contigARecord.getSequenceLength(), "Call first position greater than contig length");
        Utils.validateArg(call.getPositionB() <= contigBRecord.getSequenceLength(), "Call second position greater than contig length");
    }

    public static boolean isIntrachromosomal(final SVCallRecord call) {
        return call.getContigA().equals(call.getContigB());
    }
}
