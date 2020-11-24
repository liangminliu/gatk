package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class SVCallRecordUtils {

    public static VariantContextBuilder getVariantBuilder(final SVCallRecord call,
                                                          final Collection<String> samples) {
        Utils.nonNull(call);
        Utils.nonNull(samples);
        Utils.containsNoNull(samples, "Null sample found");
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final Allele refAllele = Allele.REF_N;
        final VariantContextBuilder builder = new VariantContextBuilder(call.getId(), call.getContig(), call.getStart(),
                call.getEnd(), Lists.newArrayList(refAllele, altAllele));
        builder.id(call.getId());
        builder.attribute(VCFConstants.END_KEY, call.getEnd());
        builder.attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, call.getContig2());
        builder.attribute(GATKSVVCFConstants.END2_ATTRIBUTE, call.getPosition2());
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        builder.attribute(GATKSVVCFConstants.SVTYPE, call.getType());
        builder.attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, getStrandString(call));
        builder.attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, call.getAlgorithms());
        final List<Genotype> genotypes = new ArrayList<>();
        for (final String sample : samples) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            if (call.getSamples().contains(sample)) {
                genotypeBuilder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
            } else {
                genotypeBuilder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
            }
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder;
    }

    public static VariantContextBuilder getVariantWithEvidenceBuilder(final SVCallRecordWithEvidence call,
                                                                      final Collection<String> samples) {
        final VariantContextBuilder builder = getVariantBuilder(call, samples);
        final SplitReadSite startSplitReadCounts = getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getStart());
        final SplitReadSite endSplitReadCounts = getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getEnd());
        final Map<String,Integer> discordantPairCounts = getDiscordantPairCountsMap(call.getDiscordantPairs());
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size());
        for (final Genotype genotype : genotypes) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            final int startCount = startSplitReadCounts.hasSample(sample) ? startSplitReadCounts.getCount(sample) : 0;
            final int endCount = endSplitReadCounts.hasSample(sample) ? startSplitReadCounts.getCount(sample) : 0;
            genotypeBuilder.attribute(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, startCount);
            genotypeBuilder.attribute(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, endCount);
            genotypeBuilder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, discordantPairCounts.getOrDefault(sample, 0));
            if (call.getSamples().contains(sample)) {
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
        return getStrandString(call.getStrand1()) + getStrandString(call.getStrand2());
    }

    private static String getStrandString(final boolean strand) {
        return strand ? SVCallRecordCodec.STRAND_PLUS : SVCallRecordCodec.STRAND_MINUS;
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

    public static int compareRecordsBySiteInfo(final SVCallRecord first, final SVCallRecord second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        final Comparator<SimpleInterval> intervalComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        // First breakpoint
        final int compareStart = intervalComparator.compare(first.getPosition1AsInterval(), second.getPosition1AsInterval());
        if (compareStart != 0) return compareStart;
        final int compareStartStrand = Boolean.compare(first.getStrand1(), first.getStrand1());
        if (compareStartStrand != 0) return compareStartStrand;
        // Second breakpoint
        final int compareEnd = intervalComparator.compare(first.getPosition2AsInterval(), second.getPosition2AsInterval());
        if (compareEnd != 0) return compareEnd;
        final int compareEndStrand = Boolean.compare(first.getStrand2(), first.getStrand2());
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

    public static Stream<SVCallRecord> convertInversionsToBreakends(final SVCallRecord call) {
        if (!call.getType().equals(StructuralVariantType.INV)) {
            return Stream.of(call);
        }
        final SVCallRecord positiveBreakend = new SVCallRecord(call.getId(), call.getContig(),
                call.getStart(), call.getEnd(), true, true, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes());
        final SVCallRecord negativeBreakend = new SVCallRecord(call.getId(), call.getContig(),
                call.getStart(), call.getEnd(), false, false, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes());
        return Stream.of(positiveBreakend, negativeBreakend);
    }

    public static <T extends Locatable> List<T> deduplicateLocatables(final List<T> items,
                                                                      final SAMSequenceDictionary dictionary,
                                                                      final BiPredicate<T,T> itemsAreIdenticalFunction,
                                                                      final Function<Collection<T>,T> deduplicateIdenticalItemsFunction) {
        final List<T> sortedItems = IntervalUtils.sortLocatablesBySequenceDictionary(items, dictionary);
        final List<T> deduplicatedList = new ArrayList<>();
        int i = 0;
        while (i < sortedItems.size()) {
            final T record = sortedItems.get(i);
            int j = i + 1;
            final Collection<Integer> identicalItemIndexes = new ArrayList<>();
            while (j < sortedItems.size() && record.getStart() == sortedItems.get(j).getStart()) {
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
}
