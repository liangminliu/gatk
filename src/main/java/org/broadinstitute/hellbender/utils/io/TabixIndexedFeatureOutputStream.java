package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.PositionalOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.util.function.Function;

public final class TabixIndexedFeatureOutputStream<F extends Feature> implements FeatureOutputStream<F> {

    private static final String NEWLINE_CHARACTER = "\n";

    private final PositionalOutputStream outputStream;
    private final IndexCreator indexCreator;
    private final Function<F, String> encoder;

    public TabixIndexedFeatureOutputStream(final GATKPath file, final FeatureCodec codec,
                                           final Function<F, String> encoder, final SAMSequenceDictionary dictionary,
                                           final int compressionLevel) {
        Utils.nonNull(file);
        Utils.nonNull(codec);
        Utils.nonNull(encoder);
        Utils.nonNull(dictionary);
        Utils.validateArg(IOUtil.hasBlockCompressedExtension(file.toPath()), "File must must have a gzip extension");
        outputStream = new PositionalOutputStream(new BlockCompressedOutputStream(file.toString(), compressionLevel));
        indexCreator = new TabixIndexCreator(dictionary, codec.getTabixFormat());
        this.encoder = encoder;
    }

    public void writeHeader(final String header) {
        Utils.nonNull(header);
        try {
            outputStream.write(header.getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing header", e);
        }
    }

    public void add(final F feature) {
        if (indexCreator != null) {
            indexCreator.addFeature(feature, outputStream.getPosition());
        }
        try {
            outputStream.write((encoder.apply(feature) + NEWLINE_CHARACTER).getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing record", e);
        }
    }

    @Override
    public void close() {
        try {
            outputStream.close();
        } catch (final IOException e) {
            throw new GATKException("Error closing output", e);
        }
    }
}
