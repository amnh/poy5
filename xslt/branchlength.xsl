<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    version="1.0">

<!-- All branches of equal length -->
<xsl:template name="equal">
    <xsl:param name="node"/>
    <xsl:param name="depth"/>
    <xsl:value-of select="$depth + 1"/>
</xsl:template>

<!-- All branches of different length -->
<xsl:template name="length">
    <xsl:param name="node"/>
    <xsl:param name="depth"/>
    <xsl:variable name="newdepth">
        <xsl:call-template name="sum_characters">
            <xsl:with-param name="nodes" select="$node/Sequence"/>
        </xsl:call-template>
    </xsl:variable>
    <xsl:value-of select="$depth + (0.5 * $newdepth)"/>
</xsl:template>

<!-- A function to calculate the lenght of a branch. -->
<xsl:template name="sum_characters">
    <xsl:param name="nodes"/>
    <xsl:choose>
        <xsl:when test="$nodes">
            <xsl:variable name="rest">
                <xsl:call-template name="sum_characters">
                    <xsl:with-param name="nodes" select="$nodes[position () > 1]"/>
                </xsl:call-template>
            </xsl:variable>
            <xsl:variable name="rest1">
                <xsl:choose>
                    <xsl:when test="starts-with($nodes[1]/@Class,'Single')">
                        <xsl:value-of select="substring-after ($nodes[1]/@Cost, '- ')"/>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:value-of select="'0'"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:variable>
            <xsl:value-of select="$rest + $rest1"/>
        </xsl:when>
        <xsl:otherwise>
            <xsl:value-of select="'1'"/>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

</xsl:transform>
