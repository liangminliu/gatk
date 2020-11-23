package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVCallRecordCodec extends AsciiFeatureCodec<SVCallRecord> {

    public static final String FORMAT_SUFFIX = ".sv.tsv.gz";
    public static final String COL_DELIMITER = "\t";
    public static String STRAND_PLUS = "+";
    public static String STRAND_MINUS = "-";
    public static final TabixFormat SV_CALL_RECORD_FORMAT = new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 3, '#', 0);
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    public SVCallRecordCodec() {
        super(SVCallRecord.class);
    }

    @Override
    public SVCallRecord decode(final String line) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 11) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        return new SVCallRecord(
                null,
                tokens.get(0),
                Integer.parseUnsignedInt(tokens.get(1)) + 1, // Convert to 1-based indexing
                Integer.parseUnsignedInt(tokens.get(2)) + 1,
                tokens.get(3).equals(STRAND_PLUS),
                tokens.get(4),
                Integer.parseUnsignedInt(tokens.get(5)) + 1,
                tokens.get(6).equals(STRAND_PLUS),
                StructuralVariantType.valueOf(tokens.get(7)),
                Integer.parseInt(tokens.get(8)),
                Arrays.asList(tokens.get(9)),
                //TODO: these aren't going to load because they're "uncalled"
                Collections.singletonList(new GenotypeBuilder().name(tokens.get(10)).make())
        );
    }

    @Override
    public TabixFormat getTabixFormat() { return SV_CALL_RECORD_FORMAT; }

    @Override
    public boolean canDecode(final String path) {
        return path.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }

    public String encode(final SVCallRecord record) {
        final List<String> data = Arrays.asList(
                record.getContig(),
                Integer.toString(record.getStart() - 1), // Convert to 0-based indexing
                Integer.toString(record.getEnd() - 1), // Convert to 0-based indexing
                record.getStrand1() ? STRAND_PLUS : STRAND_MINUS,
                record.getContig2(),
                Integer.toString(record.getPosition2() - 1),
                record.getStrand2() ? STRAND_PLUS : STRAND_MINUS,
                record.getType().name(),
                Integer.toString(record.getLength()),
                String.join(",", record.getAlgorithms()),
                String.join(",", record.getSamples())
        );
        return String.join(COL_DELIMITER, data);
    }
}
