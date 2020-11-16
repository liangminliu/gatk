package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.io.PrintStream;

/**
 * Prints SV evidence records. Can be used with -L to retrieve records on a set of intervals.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Evidence file URI
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Evidence file (local)
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk PrintSVEvidence \
 *       --evidence-file gs://my-bucket/batch_name.SR.txt.gz \
 *       -L intervals.bed \
 *       -O local.SR.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Prints SV evidence records",
        oneLineSummary = "Prints SV evidence records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class PrintSVEvidence extends FeatureWalker<Feature> {

    public static final String EVIDENCE_FILE_NAME = "evidence-file";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";
    public static final String SKIP_HEADER_NAME = "skip-header";

    @Argument(
            doc = "Input file URI with extension '.SR.txt', '.PE.txt', '.BAF.txt', or '.RD.txt' (may be gzipped).",
            fullName = EVIDENCE_FILE_NAME
    )
    private GATKPath inputFilePath;

    @Argument(
            doc = "Output file. Filenames ending in '.gz' will be block compressed.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME
    )
    private int compressionLevel = 5;

    @Argument(
            doc = "Skip printing the header",
            fullName = SKIP_HEADER_NAME
    )
    private boolean skipHeader = false;

    private PrintStream printStream;
    private IndexCreator indexCreator;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(BafEvidence.class) || featureType.equals(DepthEvidence.class)
                || featureType.equals(DiscordantPairEvidence.class) || featureType.equals(SplitReadEvidence.class);
    }

    @Override
    public GATKPath getDrivingFeaturesPath() {
        return inputFilePath;
    }

    @Override
    public void onTraversalStart() {
        if (IOUtil.hasBlockCompressedExtension(outputFile.toPath())) {
            printStream = new PrintStream(new BlockCompressedOutputStream(outputFile.toString(), compressionLevel));
            indexCreator = new TabixIndexCreator(getBestAvailableSequenceDictionary(), )
        } else {
            try {
                printStream = new PrintStream(outputFile.toString());
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not create output file", e);
            }
        }
        if (!skipHeader) {
            doHeader();
        }
    }

    private void doHeader() {
        final Object header = getDrivingFeaturesHeader();
        if (header != null) {
            if (header instanceof String) {
                printStream.println((String) header);
            } else {
                throw new GATKException.ShouldNeverReachHereException("Expected header object of type " + String.class.getSimpleName());
            }
        } else {
            logger.info("Header not found");
        }
    }

    @Override
    public void apply(final Feature feature,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        printStream.println(feature.toString());
    }

    @Override
    public Object onTraversalSuccess() {
        printStream.close();
        return null;
    }
}
