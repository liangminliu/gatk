package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.PrintStream;
import java.util.function.Function;

public final class UncompressedFeatureOutputStream<F extends Feature> implements FeatureOutputStream<F> {

    private final PrintStream outputStream;

    public UncompressedFeatureOutputStream(final GATKPath file, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(file);
        Utils.nonNull(dictionary);
        Utils.validateArg(IOUtil.hasBlockCompressedExtension(file.toPath()), "File must must have a gzip extension");
        outputStream = new PrintStream(file.getOutputStream());
    }

    public void writeHeader(final String header) {
        Utils.nonNull(header);
        try {
            outputStream.write(header.getBytes());
        } catch (final IOException e) {
            throw new GATKException("Error writing header", e);
        }
    }

    public void add(final F feature, final Function<F, String> encoder) {
        outputStream.println(encoder.apply(feature));
    }

    @Override
    public void close() {
        outputStream.close();
    }
}
