package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.CoordMath;

/**
 * Any class that has a 2-dimensional mapping onto the genome should implement Locatable2D
 * positions should be reported as 1-based and closed at both ends
 *
 */
public interface Locatable2D {

    /**
     * Gets the contig name for first coordinate
     * @return name of the contig
     */
    String getContigA();

    /**
     * @return 1-based first position
     */
    int getPositionA();

    /**
     * Gets the contig name for second coordinate
     * @return name of the contig
     */
    String getContigB();

    /**
     * @return 1-based second position
     */
    int getPositionB();

    /**
     * @return number of bases of reference covered by this interval, or 0 if coordinates on different contigs
     */
    default int getLengthOnReference() {
        return isIntrachromosomal() ? CoordMath.getLength(getPositionA(), getPositionB()) : 0;
    }

    /**
     * @return true if positions are on the same contig
     */
    default boolean isIntrachromosomal() {
        return (getContigA().equals(getContigB()));
    }
}
