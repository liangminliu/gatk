package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;

import java.util.*;
import java.util.stream.Collectors;

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
        builder.attribute(SVCluster.STRANDS_ATTRIBUTE, getStrandString(call));
        builder.attribute(SVCluster.ALGORITHMS_ATTRIBUTE, call.getAlgorithms());
        final List<Genotype> genotypes = new ArrayList<>();
        for (final String sample : samples) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            if (call.getSamples().contains(sample)) {
                genotypeBuilder.attribute(SVCluster.RAW_CALL_ATTRIBUTE, SVCluster.RAW_CALL_ATTRIBUTE_TRUE);
            } else {
                genotypeBuilder.attribute(SVCluster.RAW_CALL_ATTRIBUTE, SVCluster.RAW_CALL_ATTRIBUTE_FALSE);
            }
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder;
    }

    public static VariantContextBuilder getVariantWithEvidenceBuilder(final SVCallRecordWithEvidence call,
                                                                      final Collection<String> samples) {
        final VariantContextBuilder builder = getVariantBuilder(call, samples);
        final Map<String,Integer> startSplitReadCounts = getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getStart());
        final Map<String,Integer> endSplitReadCounts = getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getEnd());
        final Map<String,Integer> discordantPairCounts = getDiscordantPairCountsMap(call.getDiscordantPairs());
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size());
        for (final Genotype genotype : genotypes) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            genotypeBuilder.attribute(SVCluster.START_SPLIT_READ_COUNT_ATTRIBUTE, startSplitReadCounts.getOrDefault(sample, 0));
            genotypeBuilder.attribute(SVCluster.END_SPLIT_READ_COUNT_ATTRIBUTE, endSplitReadCounts.getOrDefault(sample, 0));
            genotypeBuilder.attribute(SVCluster.DISCORDANT_PAIR_COUNT_ATTRIBUTE, discordantPairCounts.getOrDefault(sample, 0));
            if (call.getSamples().contains(sample)) {
                genotypeBuilder.attribute(SVCluster.RAW_CALL_ATTRIBUTE, SVCluster.RAW_CALL_ATTRIBUTE_TRUE);
            } else {
                genotypeBuilder.attribute(SVCluster.RAW_CALL_ATTRIBUTE, SVCluster.RAW_CALL_ATTRIBUTE_FALSE);
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

    private static Map<String,Integer> getSplitReadCountsAtPosition(final List<SplitReadSite> sites, final int pos) {
        Utils.nonNull(sites);
        Utils.validateArg(pos > 0, "Non-positive position");
        if (sites.stream().map(SplitReadSite::getPosition).distinct().count() != sites.size()) {
            throw new IllegalArgumentException("Sites did not have unique positions");
        }
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .map(SplitReadSite::getSampleCountsMap)
                .findAny()
                .orElse(Collections.emptyMap());
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
}
