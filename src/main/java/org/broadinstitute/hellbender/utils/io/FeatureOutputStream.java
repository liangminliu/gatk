package org.broadinstitute.hellbender.utils.io;

import htsjdk.tribble.Feature;

import java.io.Closeable;
import java.util.function.Function;

public interface FeatureOutputStream<F extends Feature> extends Closeable {
    void writeHeader(final String header);
    void add(final F feature, final Function<F, String> encoder);
    void close();
}
