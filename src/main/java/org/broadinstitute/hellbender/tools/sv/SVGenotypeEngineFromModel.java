package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

public class SVGenotypeEngineFromModel extends SVGenotypeEngine {

    private static final String COLUMN_SEPARATOR = "\t";
    private static final String FIRST_DIM_SEPARATOR = ";";
    private static final String SECOND_DIM_SEPARATOR = ",";

    public static int getNeutralCopyNumber(final Genotype genotype) {
        Utils.validateArg(genotype.hasExtendedAttribute(SVGenotypeEngine.NEUTRAL_COPY_NUMBER_KEY),
                "Genotype missing format field " + SVGenotypeEngine.NEUTRAL_COPY_NUMBER_KEY
                        + " for sample " + genotype.getSampleName());
        Utils.nonNull(genotype);
        return Integer.valueOf(genotype.getExtendedAttribute(SVGenotypeEngine.NEUTRAL_COPY_NUMBER_KEY).toString());
    }

    public static boolean isDepthOnlyVariant(final VariantContext variant) {
        if (!variant.hasAttribute(SVCluster.ALGORITHMS_ATTRIBUTE)) {
            throw new GATKException("Variant record is missing attribute: " + SVCluster.ALGORITHMS_ATTRIBUTE);
        }
        final List<String> algorithms = variant.getAttributeAsStringList(SVCluster.ALGORITHMS_ATTRIBUTE, "");
        return algorithms.size() == 1 && algorithms.contains(SVCluster.DEPTH_ALGORITHM);
    }

    public static List<VCFHeaderLine> getVcfHeaderMetadata() {
        final List<VCFHeaderLine> lines = new ArrayList<>();
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY, true));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY, true));
        lines.add(new VCFFormatHeaderLine(SVGenotypeEngine.COPY_NUMBER_FIELD, 1, VCFHeaderLineType.Integer, "Copy number"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.PAIRED_END_PROB_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read support probability"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.FIRST_SPLIT_READ_PROB_FIELD, 1, VCFHeaderLineType.Float, "First split read support probability"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.SECOND_SPLIT_READ_PROB_FIELD, 1, VCFHeaderLineType.Float, "Second read support probability"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.PAIRED_END_BACKGROUND_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read mean background rate"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.FIRST_SPLIT_READ_BACKGROUND_FIELD, 1, VCFHeaderLineType.Float, "First split read mean background rate"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.SECOND_SPLIT_READ_BACKGROUND_FIELD, 1, VCFHeaderLineType.Float, "Second split read mean background rate"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.PAIRED_END_MEAN_BIAS_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read mean bias"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.FIRST_SPLIT_READ_MEAN_BIAS_FIELD, 1, VCFHeaderLineType.Float, "First split read mean bias"));
        lines.add(new VCFInfoHeaderLine(SVGenotypeEngine.SECOND_SPLIT_READ_MEAN_BIAS_FIELD, 1, VCFHeaderLineType.Float, "Second split read mean bias"));
        return lines;
    }

    public VariantContext genotypeFromModel(final VariantContext variant, final String modelOutputLine, final List<String> modelSampleList) {
        Utils.nonNull(variant);
        Utils.nonNull(modelOutputLine);

        final VariantOutput modelOutput = parseVariantOutput(modelOutputLine, modelSampleList.size());
        if (!modelOutput.getId().equals(variant.getID())) {
            throw new UserException.BadInput("Model and VCF record IDs did not match: " + modelOutput.getId() + ", " + variant.getID());
        }

        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        builder.attribute(SVGenotypeEngine.PAIRED_END_PROB_FIELD, modelOutput.getP_m_pe());
        builder.attribute(SVGenotypeEngine.FIRST_SPLIT_READ_PROB_FIELD, modelOutput.getP_m_sr1());
        builder.attribute(SVGenotypeEngine.SECOND_SPLIT_READ_PROB_FIELD, modelOutput.getP_m_sr2());
        builder.attribute(SVGenotypeEngine.PAIRED_END_BACKGROUND_FIELD, modelOutput.getEps_pe());
        builder.attribute(SVGenotypeEngine.FIRST_SPLIT_READ_BACKGROUND_FIELD, modelOutput.getEps_sr1());
        builder.attribute(SVGenotypeEngine.SECOND_SPLIT_READ_BACKGROUND_FIELD, modelOutput.getEps_sr2());
        builder.attribute(SVGenotypeEngine.PAIRED_END_MEAN_BIAS_FIELD, modelOutput.getPhi_pe());
        builder.attribute(SVGenotypeEngine.FIRST_SPLIT_READ_MEAN_BIAS_FIELD, modelOutput.getPhi_sr1());
        builder.attribute(SVGenotypeEngine.SECOND_SPLIT_READ_MEAN_BIAS_FIELD, modelOutput.getPhi_sr2());

        final int numSamples = modelSampleList.size();
        final StructuralVariantType svType = variant.getStructuralVariantType();
        final List<Genotype> newGenotypes = new ArrayList<>(variant.getNSamples());
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            final String sample = modelSampleList.get(sampleIndex);
            final Genotype genotype = variant.getGenotype(sample);
            if (!sample.equals(genotype.getSampleName())) {
                throw new UserException.BadInput("Model and VCF samples do not match");
            }
            final double[] genotypeProbs = modelOutput.getSampleFrequencies(sampleIndex);
            newGenotypes.add(genotypeFromGivenProbs(genotype, svType, genotypeProbs));
        }
        builder.genotypes(newGenotypes);
        final double log10ProbNoVariant = SVGenotypeEngine.calculateLog10PNoError(newGenotypes);
        builder.log10PError(log10ProbNoVariant);
        return builder.make();
    }

    public static VariantOutput parseVariantOutput(final String line, final int numSamples) {
        final String[] values = line.trim().split(COLUMN_SEPARATOR);
        final String id = values[0];
        final String[] freqStringArray = values[1].split(FIRST_DIM_SEPARATOR);
        if (freqStringArray.length != numSamples) {
            throw new UserException.BadInput("Genotype frequencies did not match sample list size");
        }
        final List<double[]> freqList = new ArrayList<>(freqStringArray.length);
        for (int i = 0; i < numSamples; i++) {
            final String[] sampleFreqStringArray = freqStringArray[i].split(SECOND_DIM_SEPARATOR);
            final double[] sampleFreq = new double[sampleFreqStringArray.length];
            for (int j = 0; j < sampleFreq.length; j++) {
                sampleFreq[j] = Double.parseDouble(sampleFreqStringArray[j]);
            }
            freqList.add(sampleFreq);
        }

        final double p_m_pe = Double.parseDouble(values[2]);
        final double p_m_sr1 = Double.parseDouble(values[3]);
        final double p_m_sr2 = Double.parseDouble(values[4]);
        final double eps_pe = Double.parseDouble(values[5]);
        final double eps_sr1 = Double.parseDouble(values[6]);
        final double eps_sr2 = Double.parseDouble(values[7]);
        final double phi_pe = Double.parseDouble(values[8]);
        final double phi_sr1 = Double.parseDouble(values[9]);
        final double phi_sr2 = Double.parseDouble(values[10]);

        return new VariantOutput(
                id,
                freqList,
                p_m_pe,
                p_m_sr1,
                p_m_sr2,
                eps_pe,
                eps_sr1,
                eps_sr2,
                phi_pe,
                phi_sr1,
                phi_sr2
        );
    }

    private static final class VariantOutput {
        private final String id;
        private final List<double[]> frequencies;
        private final double p_m_pe;
        private final double p_m_sr1;
        private final double p_m_sr2;
        private final double eps_pe;
        private final double eps_sr1;
        private final double eps_sr2;
        private final double phi_pe;
        private final double phi_sr1;
        private final double phi_sr2;

        public VariantOutput(final String id,
                             final List<double[]> frequencies,
                             final double p_m_pe,
                             final double p_m_sr1,
                             final double p_m_sr2,
                             final double eps_pe,
                             final double eps_sr1,
                             final double eps_sr2,
                             final double phi_pe,
                             final double phi_sr1,
                             final double phi_sr2) {
            this.id = id;
            this.frequencies = frequencies;
            this.p_m_pe = p_m_pe;
            this.p_m_sr1 = p_m_sr1;
            this.p_m_sr2 = p_m_sr2;
            this.eps_pe = eps_pe;
            this.eps_sr1 = eps_sr1;
            this.eps_sr2 = eps_sr2;
            this.phi_pe = phi_pe;
            this.phi_sr1 = phi_sr1;
            this.phi_sr2 = phi_sr2;
        }

        public int getNumGenotypes() {
            if (frequencies.isEmpty()) return 0;
            final int max = frequencies.stream().mapToInt(f -> f.length).max().getAsInt();
            final int min = frequencies.stream().mapToInt(f -> f.length).min().getAsInt();
            if (max != min) {
                throw new UserException.BadInput("Genotype frequency arrays not of uniform size");
            }
            return max;
        }

        public String getId() {
            return id;
        }

        public double[] getSampleFrequencies(final int sampleIndex) {
            Utils.validateArg(sampleIndex >= 0 && sampleIndex < frequencies.size(), "Invalid sample index: " + sampleIndex);
            return frequencies.get(sampleIndex);
        }

        public double getP_m_pe() {
            return p_m_pe;
        }

        public double getP_m_sr1() {
            return p_m_sr1;
        }

        public double getP_m_sr2() {
            return p_m_sr2;
        }

        public double getEps_pe() {
            return eps_pe;
        }

        public double getEps_sr1() {
            return eps_sr1;
        }

        public double getEps_sr2() {
            return eps_sr2;
        }

        public double getPhi_pe() {
            return phi_pe;
        }

        public double getPhi_sr1() {
            return phi_sr1;
        }

        public double getPhi_sr2() {
            return phi_sr2;
        }
    }
}
