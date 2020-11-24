package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;
import java.util.Objects;

public class SVCallRecordWithEvidence extends SVCallRecord {

    private final List<SplitReadSite> startSplitReadSites;
    private final List<SplitReadSite> endSplitReadSites;
    private final List<DiscordantPairEvidence> discordantPairs;
    private final CopyNumberPosteriorDistribution copyNumberDistribution;

    public SVCallRecordWithEvidence(final SVCallRecord record) {
        super(record.getId(), record.getContig(), record.getStart(), record.getEnd(), record.getStrand1(), record.getContig2(),
                record.getPosition2(), record.getStrand2(), record.getType(), record.getLength(), record.getAlgorithms(),
                record.getGenotypes());
        this.startSplitReadSites = Collections.emptyList();
        this.endSplitReadSites = Collections.emptyList();
        this.discordantPairs = Collections.emptyList();
        this.copyNumberDistribution = null;
    }

    public SVCallRecordWithEvidence(final String id,
                                    final String contig,
                                    final int start,
                                    final int end,
                                    final boolean strand1,
                                    final boolean strand2,
                                    final StructuralVariantType type,
                                    final int length,
                                    final List<String> algorithms,
                                    final List<Genotype> genotypes,
                                    final List<SplitReadSite> startSplitReadSites,
                                    final List<SplitReadSite> endSplitReadSites,
                                    final List<DiscordantPairEvidence> discordantPairs,
                                    final CopyNumberPosteriorDistribution copyNumberDistribution) {
        this(id, contig, start, end, strand1, contig, end, strand2, type, length, algorithms, genotypes,
                startSplitReadSites, endSplitReadSites, discordantPairs, copyNumberDistribution);
    }

    public SVCallRecordWithEvidence(final String id,
                                    final String startContig,
                                    final int position1,
                                    final int end1,
                                    final boolean strand1,
                                    final String contig2,
                                    final int position2,
                                    final boolean strand2,
                                    final StructuralVariantType type,
                                    final int length,
                                    final List<String> algorithms,
                                    final List<Genotype> genotypes,
                                    final List<SplitReadSite> startSplitReadSites,
                                    final List<SplitReadSite> endSplitReadSites,
                                    final List<DiscordantPairEvidence> discordantPairs,
                                    final CopyNumberPosteriorDistribution copyNumberDistribution) {
        super(id, startContig, position1, end1, strand1, contig2, position2, strand2, type, length, algorithms, genotypes);
        Utils.nonNull(startSplitReadSites);
        Utils.nonNull(endSplitReadSites);
        Utils.nonNull(discordantPairs);
        Utils.containsNoNull(startSplitReadSites, "Encountered null in start split reads");
        Utils.containsNoNull(endSplitReadSites, "Encountered null in end split reads");
        Utils.containsNoNull(discordantPairs, "Encountered null in discordant pairs");
        this.startSplitReadSites = startSplitReadSites;
        this.endSplitReadSites = endSplitReadSites;
        this.discordantPairs = discordantPairs;
        this.copyNumberDistribution = copyNumberDistribution;
    }

    public List<DiscordantPairEvidence> getDiscordantPairs() {
        return discordantPairs;
    }

    public List<SplitReadSite> getStartSplitReadSites() {
        return startSplitReadSites;
    }

    public List<SplitReadSite> getEndSplitReadSites() {
        return endSplitReadSites;
    }

    public CopyNumberPosteriorDistribution getCopyNumberDistribution() { return copyNumberDistribution; }

    @Override
    public boolean equals(final Object obj) {
        if (!super.equals(obj)) {
            return false;
        }
        if (this.getClass() != obj.getClass()) {
            return false;
        }
        final SVCallRecordWithEvidence b = (SVCallRecordWithEvidence)obj;
        boolean areEqual = this.getDiscordantPairs().containsAll(b.getDiscordantPairs());
        areEqual &= b.getDiscordantPairs().containsAll(this.getDiscordantPairs());
        areEqual &= this.getEndSplitReadSites().containsAll(b.getEndSplitReadSites());
        areEqual &= b.getEndSplitReadSites().containsAll(this.getEndSplitReadSites());
        areEqual &= this.getStartSplitReadSites().containsAll(b.getStartSplitReadSites());
        areEqual &= b.getStartSplitReadSites().containsAll(this.getStartSplitReadSites());
        areEqual &= b.getCopyNumberDistribution() == this.getCopyNumberDistribution()
                || (this.getCopyNumberDistribution() != null && this.getCopyNumberDistribution().equals(b.getCopyNumberDistribution()));
        return areEqual;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), discordantPairs, endSplitReadSites, startSplitReadSites, copyNumberDistribution);
    }

}
