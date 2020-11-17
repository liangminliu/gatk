package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.PrintStream;
import java.util.function.Function;

public final class UncompressedFeatureOutputStream<F extends Feature> implements FeatureOutputStream<F> {

    private final PrintStream outputStream;
    private final Function<F, String> encoder;

    public UncompressedFeatureOutputStream(final GATKPath file, final Function<F, String> encoder,
                                           final SAMSequenceDictionary dictionary) {
        Utils.nonNull(file);
        Utils.nonNull(encoder);
        Utils.nonNull(dictionary);
        outputStream = new PrintStream(file.getOutputStream());
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
        outputStream.println(encoder.apply(feature));
    }

    @Override
    public void close() {
        outputStream.close();
    }
}
