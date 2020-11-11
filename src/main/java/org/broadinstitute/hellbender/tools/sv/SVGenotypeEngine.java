package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class SVGenotypeEngine {

    public static final String COPY_NUMBER_LOG_POSTERIORS_KEY = "CNLP";
    public static final String NEUTRAL_COPY_NUMBER_KEY = "NCN";
    public static final String COPY_NUMBER_FIELD = "CN";
    public static final String PAIRED_END_PROB_FIELD = "PPE";
    public static final String FIRST_SPLIT_READ_PROB_FIELD = "PSR1";
    public static final String SECOND_SPLIT_READ_PROB_FIELD = "PSR2";
    public static final String PAIRED_END_BACKGROUND_FIELD = "EPE";
    public static final String FIRST_SPLIT_READ_BACKGROUND_FIELD = "ESR1";
    public static final String SECOND_SPLIT_READ_BACKGROUND_FIELD = "ESR2";
    public static final String PAIRED_END_MEAN_BIAS_FIELD = "PHI_PE";
    public static final String FIRST_SPLIT_READ_MEAN_BIAS_FIELD = "PHI_SR1";
    public static final String SECOND_SPLIT_READ_MEAN_BIAS_FIELD = "PHI_SR2";
    public static final Allele DEFAULT_REF_ALLELE = Allele.REF_N;
    public static final String BND_SYMBOLIC_ALLELE = "<" + StructuralVariantType.BND.name() + ">";
    public static final List<StructuralVariantType> CNV_TYPES = Lists.newArrayList(StructuralVariantType.DEL, StructuralVariantType.DUP);

    public static double calculateLog10PNoError(final Collection<Genotype> genotypes) {
        double log10ProbNoVariant = 0;
        for (final Genotype genotype : genotypes) {
            final int[] genotypeQuals = genotype.getPL();
            if (genotypeQuals.length > 0) {
                log10ProbNoVariant += QualityUtils.qualToErrorProbLog10(genotypeQuals[0]);
            }
        }
        return log10ProbNoVariant;
    }

    public static Genotype buildGenotypeFromQuals(final Genotype genotype, final StructuralVariantType svType,
                                          final int[] genotypeQuals, final int genotypeIndex,
                                          final int genotypeQuality) {
        final Allele altAllele = Allele.create("<" + svType.name() + ">", false);
        final Allele refAllele = DEFAULT_REF_ALLELE;
        final int neutralCopyState = SVGenotypeEngineFromModel.getNeutralCopyNumber(genotype);
        // TODO: multi-allelic sites
        final List<Allele> alleles = getAlleles(svType, genotypeIndex, neutralCopyState, refAllele, altAllele);
        final int copyNumber = getGenotypeCopyNumber(svType, genotypeIndex, neutralCopyState);

        final GenotypeBuilder builder = new GenotypeBuilder(genotype);
        builder.alleles(alleles);
        builder.attribute(COPY_NUMBER_FIELD, copyNumber);
        builder.GQ(genotypeQuality);
        builder.PL(genotypeQuals);
        return builder.make();
    }

    private static Integer getGenotypeCopyNumber(final StructuralVariantType svType, final int genotypeIndex, final int neutralCopyState) {
        if (svType.equals(StructuralVariantType.INS)
                || svType.equals(StructuralVariantType.INV)
                || svType.equals(StructuralVariantType.BND)) {
            return null;
        } else if (svType.equals(StructuralVariantType.DEL)) {
            return neutralCopyState - genotypeIndex;
        } else if (svType.equals(StructuralVariantType.DUP)) {
            return neutralCopyState + genotypeIndex;
        }
        throw new UserException.BadInput("Unsupported SVTYPE: " + svType);
    }

    private static List<Allele> getAlleles(final StructuralVariantType svType, final int genotypeIndex,
                                           final int neutralCopyState, final Allele refAllele, final Allele altAllele) {
        if (svType.equals(StructuralVariantType.DEL)
                || svType.equals(StructuralVariantType.INS)
                || svType.equals(StructuralVariantType.INV)
                || svType.equals(StructuralVariantType.BND)) {
            return getNonDupAlleles(genotypeIndex, neutralCopyState, refAllele, altAllele);
        } else if (svType.equals(StructuralVariantType.DUP)) {
            return getDupAlleles(neutralCopyState);
        }
        throw new UserException.BadInput("Unsupported SVTYPE: " + svType);
    }

    private static List<Allele> getNonDupAlleles(final int genotypeIndex, final int neutralCopyState,
                                                 final Allele refAllele, final Allele altAllele) {
        final int alleleCount = getAlleleCount(neutralCopyState);
        final int numAltAlleles = Math.min(genotypeIndex, neutralCopyState);
        final List<Allele> alleles = new ArrayList<>(neutralCopyState);
        for (int i = 0; i < alleleCount - numAltAlleles; i++) {
            alleles.add(refAllele);
        }
        for (int i = 0; i < numAltAlleles; i++) {
            alleles.add(altAllele);
        }
        return alleles;
    }

    private static List<Allele> getDupAlleles(final int neutralCopyState) {
        final int alleleCount = getAlleleCount(neutralCopyState);
        final List<Allele> alleles = new ArrayList<>(alleleCount);
        for (int i = 0; i < alleleCount; i++) {
            alleles.add(Allele.NO_CALL);
        }
        return alleles;
    }

    private static int getAlleleCount(final int neutralCopyState) {
        // Can't have mixed rows with some samples missing genotypes, so clamp to a minimum of 1
        return Math.max(1, neutralCopyState);
    }
}
