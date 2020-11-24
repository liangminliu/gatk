package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecord implements Feature {

    private final String id;
    private final String contig1;
    private final int position1;
    private final int end1;  // END / stop position (must be >= start)
    private final boolean strand1;
    private final String contig2;
    private final int position2;
    private final boolean strand2;
    private final StructuralVariantType type;
    private int length;
    private final List<String> algorithms;
    private final List<Genotype> genotypes;
    private LinkedHashSet<String> samples;

    private final static List<String> nonDepthCallerAttributes = Arrays.asList(
            VCFConstants.END_KEY,
            GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE,
            GATKSVVCFConstants.STRANDS_ATTRIBUTE,
            GATKSVVCFConstants.SVLEN,
            GATKSVVCFConstants.SVTYPE
    );

    public static SVCallRecord create(final VariantContext variant) {
        Utils.nonNull(variant);
        Utils.validate(variant.getAttributes().keySet().containsAll(nonDepthCallerAttributes), "Call is missing attributes");
        final String id = variant.getID();
        final String contig1 = variant.getContig();
        final int position1 = variant.getStart();
        final int end1 = variant.getEnd();
        final String contig2;
        final int position2;

        // If END2 and CONTIG2 are both defined, use those.
        // If neither is defined, use start contig and position.
        // If only CONTIG2 is defined, END2 is taken as END
        // Having only END2 is unacceptable
        final boolean hasContig2 = variant.hasAttribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
        final boolean hasEnd2 = variant.hasAttribute(GATKSVVCFConstants.END2_ATTRIBUTE);
        if (hasContig2 && hasEnd2) {
            contig2 = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            position2 = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, end1);
        } else if (!hasContig2 && !hasEnd2) {
            contig2 = contig1;
            position2 = position1;
        } else if (hasContig2) {
            contig2 = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            position2 = end1;
        } else {
            throw new UserException.BadInput("Attribute " + GATKSVVCFConstants.END2_ATTRIBUTE +
                    " cannot be defined without " + GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
        }

        final StructuralVariantType type = variant.getStructuralVariantType();
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE), "Attribute " + GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE + " is required");
        final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE), "Attribute " + GATKSVVCFConstants.STRANDS_ATTRIBUTE + " is required");
        final String strands = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, null);
        if (strands.length() != 2) {
            throw new IllegalArgumentException("Strands field is not 2 characters long");
        }
        final String strand1Char = strands.substring(0, 1);
        if (!strand1Char.equals(SVCallRecordCodec.STRAND_PLUS) && !strand1Char.equals(SVCallRecordCodec.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final String strand2Char = strands.substring(1, 2);
        if (!strand2Char.equals(SVCallRecordCodec.STRAND_PLUS) && !strand2Char.equals(SVCallRecordCodec.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        final boolean strand1 = strand1Char.equals(SVCallRecordCodec.STRAND_PLUS);
        final boolean strand2 = strand2Char.equals(SVCallRecordCodec.STRAND_PLUS);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.SVLEN), "Attribute " + GATKSVVCFConstants.SVLEN + " is required");
        final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        return new SVCallRecord(id, contig1, position1, end1, strand1, contig2, position2, strand2, type, length, algorithms, variant.getGenotypes());
    }

    /**
     *
     * @param variant single-sample variant from a gCNV segments VCF
     * @param minQuality drop events with quality lower than this
     * @return
     */
    public static SVCallRecord createDepthOnlyFromGCNV(final VariantContext variant, final double minQuality) {
        Utils.nonNull(variant);

        if (variant.getGenotypes().size() == 1) {
            //only cluster good variants
            final Genotype g = variant.getGenotypes().get(0);
            if (g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))
                    || Integer.valueOf((String) g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS)) < minQuality) {
                return null;
            }
        }


        final List<String> algorithms = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);

        boolean isDel = false;
        for (final Genotype g : variant.getGenotypes()) {
            if (g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))) {
                continue;
            }
            if (variant.getReference().equals(Allele.REF_N)) {  //old segments VCFs had ref Ns and genotypes that didn't reflect ploidy accurately
                if (g.getAlleles().stream().anyMatch(a -> a.equals(GATKSVVCFConstants.DEL_ALLELE))) {
                    isDel = true;
                } else if (g.getAlleles().stream().anyMatch(a -> a.equals(GATKSVVCFConstants.DUP_ALLELE))) {
                    isDel = false;
                } else if (g.getAlleles().stream().allMatch(a -> a.isNoCall())) {
                    if (g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                        isDel = (Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) < g.getPloidy());
                    } else {
                        throw new IllegalStateException("Genotype for sample " + g.getSampleName() + " at " + variant.getContig() + ":" + variant.getStart() + " had no CN attribute and will be dropped.");
                    }
                } else {
                    throw new IllegalArgumentException("Segment VCF schema expects <DEL>, <DUP>, and no-call allele, but found "
                            + g.getAllele(0) + " at " + variant.getContig() + ":" + variant.getStart());
                }
            } else {  //for spec-compliant VCFs (i.e. with non-N ref allele) we can just trust the ALT
                isDel = (variant.getAlternateAlleles().contains(GATKSVVCFConstants.DEL_ALLELE)
                        && !variant.getAlternateAlleles().contains(GATKSVVCFConstants.DUP_ALLELE));
            }
        }

        final boolean startStrand = isDel ? true : false;
        final boolean endStrand = isDel ? false : true;
        final StructuralVariantType type;
        if (!variant.getReference().equals(Allele.REF_N) && variant.getAlternateAlleles().contains(GATKSVVCFConstants.DUP_ALLELE)
                && variant.getAlternateAlleles().contains(GATKSVVCFConstants.DEL_ALLELE)) {
            type = StructuralVariantType.CNV;
        } else {
            type = isDel ? StructuralVariantType.DEL : StructuralVariantType.DUP;
        }

        final String id = variant.getID();
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final int end = variant.getEnd();
        final int length = end - start;
        return new SVCallRecord(id, startContig, start, end, startStrand, startContig, end, endStrand, type, length, algorithms,
                new ArrayList<>(variant.getGenotypes()));
    }

    public SVCallRecord(final String id,
                        final String contig,
                        final int start,
                        final int end,
                        final boolean strand1,
                        final boolean strand2,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final List<Genotype> genotypes) {
        this(id, contig, start, end, strand1, contig, end, strand2, type, length, algorithms, genotypes);
    }

    public SVCallRecord(final String id,
                        final String contig1,
                        final int position1,
                        final int end1,
                        final boolean strand1,
                        final String contig2,
                        final int position2,
                        final boolean strand2,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final List<Genotype> genotypes) {
        Utils.nonNull(id);
        Utils.nonNull(contig1);
        Utils.nonNull(contig2);
        Utils.nonNull(type);
        Utils.nonNull(algorithms);
        Utils.nonNull(genotypes);
        Utils.nonEmpty(algorithms);
        Utils.nonEmpty(genotypes);
        Utils.containsNoNull(algorithms, "Encountered null algorithm");
        Utils.containsNoNull(genotypes, "Encountered null genotype");
        this.id = id;
        this.contig1 = contig1;
        this.position1 = position1;
        this.end1 = end1;
        this.strand1 = strand1;
        this.contig2 = contig2;
        this.position2 = position2;
        this.strand2 = strand2;
        this.type = type;
        this.length = length;
        this.algorithms = Collections.unmodifiableList(algorithms);
        this.genotypes = Collections.unmodifiableList(genotypes);
        this.samples = genotypes.stream()
                .filter(g -> !g.getType().equals(GenotypeType.NO_CALL) || g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) //skip no-calls that don't have copy number (i.e. not dupe no-calls)
                .map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    public String getId() {
        return id;
    }

    @Override
    public String getContig() {
        return contig1;
    }

    @Override
    public int getStart() {
        return position1;
    }

    public boolean getStrand1() {
        return strand1;
    }

    public String getContig2() {
        return contig2;
    }

    @Override
    public int getEnd() {
        return end1;
    }

    public int getPosition2() {
        return position2;
    }

    public boolean getStrand2() {
        return strand2;
    }

    public StructuralVariantType getType() {
        return type;
    }

    public int getLength() {
        return length;
    }

    public List<String> getAlgorithms() {
        return algorithms;
    }

    public Set<String> getSamples() {
        return samples;
    }

    public List<Genotype> getGenotypes() {
        return genotypes;
    }

    public SimpleInterval getPosition1AsInterval() {
        return new SimpleInterval(contig1, position1, position1);
    }

    public SimpleInterval getPosition2AsInterval() {
        return new SimpleInterval(contig2, end1, end1);
    }

    String prettyPrint() {
        return getContig() + ":" + getStart() + "-" + getEnd();
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null) {
            return false;
        }
        if (this.getClass() != obj.getClass()) {
            return false;
        }
        final SVCallRecord b = (SVCallRecord) obj;

        //quick check
        if (!this.getContig().equals(b.getContig())) return false;
        if (this.getStart() != b.getStart()) return false;

        boolean areEqual = this.getStrand1() == b.getStrand1();

        areEqual &= this.getId() == b.getId();
        areEqual &= this.getContig2() == b.getContig2();
        areEqual &= this.getEnd() == b.getEnd();
        areEqual &= this.getStrand2() == b.getStrand2();

        areEqual &= this.getType() == b.getType();
        areEqual &= this.getLength() == b.getLength();

        areEqual &= this.getAlgorithms().containsAll(b.getAlgorithms());
        areEqual &= b.getAlgorithms().containsAll(this.getAlgorithms());

        areEqual &= this.getSamples().containsAll(b.getSamples());
        areEqual &= b.getSamples().containsAll(this.getSamples());

        areEqual &= this.getGenotypes().containsAll(b.getGenotypes());
        areEqual &= b.getGenotypes().containsAll(this.getGenotypes());

        return areEqual;
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, algorithms, end1, contig2, strand2, genotypes,
                length, samples, position1, contig1, strand1, type);
    }
}
