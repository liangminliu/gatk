package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecord implements Locatable2D {

    public static final String STRAND_PLUS = "+";
    public static final String STRAND_MINUS = "-";

    private final String id;
    private final String contigA;
    private final int positionA;
    private final boolean strandA;
    private final String contigB;
    private final int positionB;
    private final boolean strandB;
    private final StructuralVariantType type;
    private int length;
    private final List<String> algorithms;
    private final GenotypesContext genotypes;

    private Set<String> allSamples;
    private Set<String> calledSamples;
    private Set<String> carrierSamples;

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
        final String contigA = variant.getContig();
        final int positionA = variant.getStart();
        final int end = variant.getEnd();
        final String contigB;
        final int positionB;

        // If END2 and CONTIG2 are both defined, use those.
        // If neither is defined, use start contig and position.
        // If only CONTIG2 is defined, END2 is taken as END
        // Having only END2 is unacceptable
        final boolean hasContig2 = variant.hasAttribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
        final boolean hasEnd2 = variant.hasAttribute(GATKSVVCFConstants.END2_ATTRIBUTE);
        if (hasContig2 && hasEnd2) {
            contigB = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            positionB = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, 0);
        } else if (!hasContig2 && !hasEnd2) {
            contigB = contigA;
            positionB = positionA;
        } else if (hasContig2) {
            contigB = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            positionB = end;
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
        if (!strand1Char.equals(STRAND_PLUS) && !strand1Char.equals(STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final String strand2Char = strands.substring(1, 2);
        if (!strand2Char.equals(STRAND_PLUS) && !strand2Char.equals(STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        final boolean strand1 = strand1Char.equals(STRAND_PLUS);
        final boolean strand2 = strand2Char.equals(STRAND_PLUS);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.SVLEN), "Attribute " + GATKSVVCFConstants.SVLEN + " is required");
        final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        return new SVCallRecord(id, contigA, positionA, strand1, contigB, positionB, strand2, type, length, algorithms, variant.getGenotypes());
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
        return new SVCallRecord(id, startContig, start, startStrand, startContig, end, endStrand, type, length, algorithms, variant.getGenotypes());
    }

    public SVCallRecord(final String id,
                        final String contig,
                        final int start,
                        final int end,
                        final boolean strandA,
                        final boolean strandB,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final List<Genotype> genotypes) {
        this(id, contig, start, strandA, contig, end, strandB, type, length, algorithms, genotypes);
    }

    public SVCallRecord(final String id,
                        final String contigA,
                        final int positionA,
                        final boolean strandA,
                        final String contigB,
                        final int positionB,
                        final boolean strandB,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final List<Genotype> genotypes) {
        Utils.nonNull(id);
        Utils.nonNull(contigA);
        Utils.nonNull(contigB);
        Utils.nonNull(type);
        Utils.nonNull(algorithms);
        Utils.nonNull(genotypes);
        Utils.nonEmpty(algorithms);
        Utils.nonEmpty(genotypes);
        Utils.containsNoNull(algorithms, "Encountered null algorithm");
        Utils.containsNoNull(genotypes, "Encountered null genotype");
        this.id = id;
        this.contigA = contigA;
        this.positionA = positionA;
        this.strandA = strandA;
        this.contigB = contigB;
        this.positionB = positionB;
        this.strandB = strandB;
        this.type = type;
        this.length = length;
        this.algorithms = Collections.unmodifiableList(algorithms);
        this.genotypes = GenotypesContext.copy(genotypes);
    }

    public String getId() {
        return id;
    }

    @Override
    public String getContigA() {
        return contigA;
    }

    @Override
    public int getPositionA() {
        return positionA;
    }

    @Override
    public String getContigB() {
        return contigB;
    }

    @Override
    public int getPositionB() {
        return positionB;
    }

    public boolean getStrandA() {
        return strandA;
    }

    public boolean getStrandB() {
        return strandB;
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

    public Set<String> getAllSamples() {
        if (allSamples == null) {
            allSamples = genotypes.stream().map(Genotype::getSampleName)
                    .collect(Collectors.toCollection(LinkedHashSet::new));
        }
        return allSamples;
    }

    public Set<String> getCalledSamples() {
        if (calledSamples == null) {
            calledSamples = genotypes.stream()
                    .filter(g -> !g.getType().equals(GenotypeType.NO_CALL) || g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) //skip no-calls that don't have copy number (i.e. not dupe no-calls)
                    .map(Genotype::getSampleName)
                    .collect(Collectors.toCollection(LinkedHashSet::new));
        }
        return calledSamples;
    }

    public Set<String> getCarrierSamples() {
        if (carrierSamples == null) {
            carrierSamples = genotypes.stream()
                    .filter(SVCallRecord::isCarrier)
                    .map(Genotype::getSampleName)
                    .collect(Collectors.toCollection(LinkedHashSet::new));
        }
        return carrierSamples;
    }

    public static boolean isCarrier(final Genotype g) {
        Utils.nonNull(g);
        return g.getType().equals(GenotypeType.HET) || g.getType().equals(GenotypeType.HOM_VAR);
    }

    public static boolean isRawCall(final Genotype g) {
        Utils.nonNull(g);
        return VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE) == GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE;
    }

    public List<Genotype> getGenotypes() {
        return genotypes;
    }

    public SimpleInterval getPositionAAsInterval() {
        return new SimpleInterval(contigA, positionA, positionA);
    }

    public SimpleInterval getPositionBAsInterval() {
        return new SimpleInterval(contigB, positionB, positionB);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SVCallRecord)) return false;
        SVCallRecord that = (SVCallRecord) o;
        return positionA == that.positionA &&
                strandA == that.strandA &&
                positionB == that.positionB &&
                strandB == that.strandB &&
                length == that.length &&
                id.equals(that.id) &&
                contigA.equals(that.contigA) &&
                contigB.equals(that.contigB) &&
                type == that.type &&
                algorithms.equals(that.algorithms) &&
                genotypes.equals(that.genotypes);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, contigA, positionA, strandA, contigB, positionB, strandB, type, length, algorithms, genotypes);
    }
}
