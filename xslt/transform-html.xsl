<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="xml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" indent="yes"/>

<xsl:strip-space elements="*"/>

<xsl:key name="alphabet" match="Element" use="@Value" />

    <xsl:template match="/">
        <html><head>
        <title>POY4 HTML Output</title></head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
        <link rel="stylesheet" type="text/css" media="screen" href="style.css" /><body>
        <div id="header"><span>POY</span></div>
        <xsl:apply-templates select="Diagnosis/Forest/Tree"/>
        </body></html>
    </xsl:template>

    <xsl:template match="Diagnosis/Forest/Tree">
        <div class="forest">
            <xsl:for-each select=".//Node">
                <div class="node">
                    <h1>node: <xsl:value-of select="normalize-space(@Name)"/></h1>
                    <xsl:if test="Non_Additive_Character != '' and position() != 1">
                        <xsl:for-each select="Non_Additive_Character">
                            <xsl:variable name="characterName" select="@Name"/>
                            <xsl:variable name="entityPosition" select="position()"/>
                            <xsl:variable name="ancestorValue" select="normalize-space(../../../Node/*[$entityPosition])"/>
                            <xsl:variable name="descendantValue" select="normalize-space(.)"/>
                            <xsl:if test="( $descendantValue != $ancestorValue ) and @Class='Final'">
                                <xsl:choose>
                                    <xsl:when test="translate($descendantValue,'ACGT- ','') = '' and (count(child::*) = string-length(translate($descendantValue,' ','')))">
                                        <h2><span>transformation</span></h2>
                                        <xsl:variable name="ancA" select="contains($ancestorValue, 'A')"/>
                                        <xsl:variable name="ancC" select="contains($ancestorValue, 'C')"/>
                                        <xsl:variable name="ancG" select="contains($ancestorValue, 'G')"/>
                                        <xsl:variable name="ancT" select="contains($ancestorValue, 'T')"/>
                                        <xsl:variable name="anc-" select="contains($ancestorValue, '-')"/>
                                        <xsl:variable name="descA" select="contains($descendantValue, 'A')"/>
                                        <xsl:variable name="descC" select="contains($descendantValue, 'C')"/>
                                        <xsl:variable name="descG" select="contains($descendantValue, 'G')"/>
                                        <xsl:variable name="descT" select="contains($descendantValue, 'T')"/>
                                        <xsl:variable name="desc-" select="contains($descendantValue, '-')"/>
                                        <xsl:variable name="ancSize">
                                            <xsl:choose>
                                                <xsl:when test="($ancA and $ancC) and ($ancG and $ancT)">
                                                    <xsl:text>0</xsl:text>
                                                </xsl:when>
                                                <xsl:otherwise>
                                                    <xsl:value-of select="number($ancA)+number($ancC)+number($ancG)+number($ancT)"/>
                                                </xsl:otherwise>
                                            </xsl:choose>
                                        </xsl:variable>
                                        <xsl:variable name="descSize">
                                            <xsl:choose>
                                                <xsl:when test="($descA and $descC) and ($descG and $descT)">
                                                    <xsl:text>0</xsl:text>
                                                </xsl:when>
                                                <xsl:otherwise>
                                                    <xsl:value-of select="number($descA)+number($descC)+number($descG)+number($descT)"/>
                                                </xsl:otherwise>
                                            </xsl:choose>
                                        </xsl:variable>
                                        <xsl:variable name="intersectionSize" select="
                                            number($ancA and $descA) +
                                            number($ancC and $descC) + 
                                            number($ancG and $descG) +
                                            number($ancT and $descT)"/>
                                        <table>
                                            <tr>
                                                <th>character</th>
                                                <td><xsl:value-of select="substring-before(@Name,':')"/></td>
                                            </tr><tr>
                                                <th>position</th>
                                                <td><xsl:value-of select="substring-after(@Name,':')"/></td>
                                            </tr><tr>
                                                <th>ancestor state</th>
                                                <td>
                                                    <xsl:if test="$ancA">A</xsl:if>
                                                    <xsl:if test="$ancC">C</xsl:if>
                                                    <xsl:if test="$ancG">G</xsl:if>
                                                    <xsl:if test="$ancT">T</xsl:if>
                                                    <xsl:if test="$anc-">-</xsl:if>
                                                </td>
                                            </tr><tr>
                                                <th>descendant state</th>
                                                <td>
                                                    <xsl:if test="$descA">A</xsl:if>
                                                    <xsl:if test="$descC">C</xsl:if>
                                                    <xsl:if test="$descG">G</xsl:if>
                                                    <xsl:if test="$descT">T</xsl:if>
                                                    <xsl:if test="$desc-">-</xsl:if>
                                                </td>
                                            </tr><tr>
                                                <th>type</th>
                                                <td>
                                                    <xsl:choose>
                                                        <xsl:when test="$ancSize = 0">insertion</xsl:when>
                                                        <xsl:when test="$descSize = 0">deletion</xsl:when>
                                                        <xsl:otherwise>
                                                            <xsl:choose>
                                                                <!-- when intersection is empty set or there is a state of length 3, it's not definite -->
                                                                <xsl:when test="$intersectionSize != 0 or ($ancSize = 3 or $descSize=3)">ambiguous</xsl:when>
                                                                <!-- changes between states of set size 1 -->
                                                                <xsl:when test="$ancSize = 1 and $descSize = 1">
                                                                    <xsl:choose>
                                                                        <xsl:when test="(($ancA or $ancG) and ($descC or $descT)) or
                                                                        (($ancC or $ancT) and ($descA or $descG))">transversion</xsl:when>
                                                                        <xsl:otherwise>transition</xsl:otherwise>
                                                                    </xsl:choose>
                                                                </xsl:when>
                                                                <!-- changes between states of set size 2 -->
                                                                <xsl:when test="$ancSize = 2 and $descSize = 2">
                                                                    <xsl:choose>
                                                                        <xsl:when test="($ancA and $ancG) or ($ancC and $ancT)">transversion</xsl:when>
                                                                        <xsl:otherwise>ambiguous</xsl:otherwise>
                                                                    </xsl:choose>
                                                                </xsl:when>
                                                                <!-- cases left: set sizes of 1 and 2 -->
                                                                <xsl:otherwise>
                                                                    <xsl:choose>
                                                                        <xsl:when test="
                                                                        ($ancA and not($descG)) or
                                                                        ($ancC and not($descT)) or
                                                                        ($ancG and not($descA)) or
                                                                        ($ancT and not($descC))
                                                                        ">transversion</xsl:when>
                                                                        <xsl:otherwise>ambiguous</xsl:otherwise>
                                                                    </xsl:choose>
                                                                </xsl:otherwise>
                                                            </xsl:choose>
                                                        </xsl:otherwise>
                                                    </xsl:choose>
                                                </td>
                                            </tr><tr>
                                                <th>cost</th>
                                                <td><xsl:value-of select="@Cost"/></td>
                                            </tr><tr>
                                                <th>definite</th>
                                                <td>
                                                    <xsl:choose>
                                                        <xsl:when test="$intersectionSize != 0">false</xsl:when>
                                                        <xsl:otherwise>true</xsl:otherwise>
                                                    </xsl:choose>
                                                </td>
                                            </tr>
                                        </table>
                                    </xsl:when>
                                    <xsl:otherwise>
                                        <h2><span>transformation</span></h2>
                                        <table>
                                            <tr>
                                                <th>character</th>
                                                <td><xsl:value-of select="translate(@Name,'_',' ')"/></td>
                                            </tr><tr>
                                                <th>ancestor state</th>
                                                <td><xsl:value-of select="translate($ancestorValue,'_',' ')"/></td>
                                            </tr><tr>
                                                <th>descendant state</th>
                                                <td><xsl:value-of select="translate($descendantValue,'_',' ')"/></td>
                                            </tr><tr>
                                                <th>type</th>
                                                <td>non-additive character</td>
                                            </tr><tr>
                                                <th>cost</th>
                                                <td><xsl:value-of select="@Cost"/></td>
                                            </tr><tr>
                                                <th>definite</th>
                                                <td><xsl:value-of select="@Definite"/></td>
                                            </tr>
                                        </table>
                                    </xsl:otherwise>
                                </xsl:choose>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="Additive_Character != '' and position() != 1">
                        <xsl:for-each select="Additive_Character">
                            <xsl:variable name="characterName" select="@Name"/>
                            <xsl:variable name="entityPosition" select="position()"/>
                            <xsl:variable name="ancestorMin" select="normalize-space(../../../Node/*[$entityPosition]/Min)"/>
                            <xsl:variable name="ancestorMax" select="normalize-space(../../../Node/*[$entityPosition]/Max)"/>
                            <xsl:if test="(( normalize-space(Min) != $ancestorMin ) or ( normalize-space(Max) != $ancestorMax )) and @Class='Final'">
                                <h2><span>transformation</span></h2>
                                <table>
                                    <tr>
                                        <th>character</th>
                                        <td><xsl:value-of select="translate(@Name,'_',' ')"/></td>
                                    </tr><tr>
                                        <th>ancestor state</th>
                                        <td>
                                            <xsl:choose>
                                                <xsl:when test="$ancestorMin = $ancestorMax">
                                                    <xsl:value-of select="translate($ancestorMin,'_',' ')"/>
                                                </xsl:when>
                                                <xsl:otherwise>
                                                    <xsl:value-of select="translate($ancestorMin,'_',' ')"/> - <xsl:value-of select="translate($ancestorMax,'_',' ')"/>
                                                </xsl:otherwise>
                                            </xsl:choose>
                                    </td>
                                    </tr><tr>
                                        <th>descendant state</th>
                                        <td>
                                            <xsl:choose>
                                                <xsl:when test="normalize-space(Min) = normalize-space(Max)">
                                                    <xsl:value-of select="translate(normalize-space(Min),'_',' ')"/>
                                                </xsl:when>
                                                <xsl:otherwise>
                                                    <xsl:value-of select="translate(normalize-space(Min),'_',' ')"/> - <xsl:value-of select="translate(normalize-space(Max),'_',' ')"/>
                                                </xsl:otherwise>
                                            </xsl:choose>
                                        </td>
                                    </tr><tr>
                                        <th>type</th>
                                        <td>additive character</td>
                                    </tr><tr>
                                        <th>cost</th>
                                        <td><xsl:value-of select="@Cost"/></td>
                                    </tr><tr>
                                        <th>definite</th>
                                        <td><xsl:value-of select="@Definite"/></td>
                                    </tr>
                                </table>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="Sankoff_Character != '' and position() != 1">
                        <xsl:for-each select="Sankoff_Character">
                            <xsl:variable name="characterName" select="@Name"/>
                            <xsl:variable name="entityPosition" select="position()"/>
                            <xsl:variable name="ancestorValue" select="normalize-space(../../../Node/*[$entityPosition]/Value)"/>
                            <xsl:if test="normalize-space(Value) != $ancestorValue and @Class='Final'">
                                <h2><span>transformation</span></h2>
                                <table>
                                    <tr>
                                        <th>character</th>
                                        <td><xsl:value-of select="substring-before(@Name,':')"/></td>
                                    </tr><tr>
                                        <th>position</th>
                                        <td><xsl:value-of select="substring-after(@Name,':')"/></td>
                                    </tr><tr>
                                        <th>ancestor state</th>
                                        <td><xsl:value-of select="$ancestorValue"/></td>
                                    </tr><tr>
                                        <th>descendant state</th>
                                        <td><xsl:value-of select="normalize-space(Value)"/></td>
                                    </tr><tr>
                                        <th>type</th>
                                        <td>sankoff character</td>
                                        </tr><tr>
                                        <th>cost</th>
                                        <td><xsl:value-of select="@Cost"/></td>
                                    </tr><tr>
                                        <th>definite</th>
                                        <td><xsl:value-of select="@Definite"/></td>
                                    </tr>
                                </table>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="Sequence != '' and position() != 1">
                        <xsl:for-each select="Sequence">
                            <xsl:variable name="characterName" select="@Name"/>
                            <xsl:variable name="entityPosition" select="position()"/>
                            <xsl:variable name="ancestorValue" select="normalize-space(../../../Node/*[$entityPosition])"/>
                            <xsl:if test="normalize-space(.) != $ancestorValue and @Class='Single'">
                                <h2><span>transformation</span></h2>
                                <div class="transformations"><table>
                                    <tr>
                                        <th>character</th>
                                        <td><xsl:value-of select="substring-before($characterName,':')"/></td>
                                    </tr><tr>
                                        <th>sequence</th>
                                        <td><xsl:value-of select="substring-after($characterName,':')"/></td>
                                    </tr><tr>
                                        <th>ancestor sequence</th>
                                        <td><xsl:value-of select="$ancestorValue"/></td>
                                    </tr><tr>
                                        <th>descendant sequence</th>
                                        <td><xsl:value-of select="normalize-space(.)"/></td>
                                    </tr><tr>
                                        <th>type</th>
                                        <td>sequence</td>
                                    </tr><tr>
                                        <th>cost</th>
                                        <td><xsl:value-of select="@Cost"/></td>
                                    </tr><tr>
                                        <th>definite</th>
                                        <td><xsl:value-of select="@Definite"/></td>
                                    </tr>
                                </table></div>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="Breakinv != '' and position() != 1">
                        <xsl:for-each select="Breakinv">
                            <xsl:variable name="characterName" select="@Name"/>
                            <xsl:variable name="entityPosition" select="position()"/>
                            <xsl:variable name="ancestorValue" select="normalize-space(../../../Node/*[$entityPosition])"/>
                            <xsl:if test="(normalize-space(.) != $ancestorValue) and @Class='Single'">
                                <xsl:variable name="totalCost" select="@Cost"/>
                                <xsl:variable name="rearrangementCost" select="@Rearrangment_Cost"/>
                                <h2><span>transformation</span></h2>
                                <div class="transformations"><table>
                                    <tr>
                                        <th>character</th>
                                        <td><xsl:value-of select="substring-before($characterName,':')"/></td>
                                    </tr><tr>
                                        <th>sequence</th>
                                        <td><xsl:value-of select="substring-after($characterName,':')"/></td>
                                    </tr><tr>
                                        <th>ancestor sequence</th>
                                        <td><xsl:value-of select="$ancestorValue"/></td>
                                    </tr><tr>
                                        <th>descendant sequence</th>
                                        <td><xsl:value-of select="normalize-space(.)"/></td>
                                    </tr><tr>
                                        <th>type</th>
                                        <td>breakinv</td>
                                    </tr><tr>
                                        <th>edit cost</th>
                                        <td><xsl:value-of select="$totalCost - $rearrangementCost"/></td>
                                    </tr><tr>
                                        <th>rearrangement cost</th>
                                        <td><xsl:value-of select="$rearrangementCost"/></td>
                                    </tr><tr>
                                        <th>total cost</th>
                                        <td><xsl:value-of select="$totalCost"/></td>
                                    </tr><tr>
                                        <th>definite</th>
                                        <td><xsl:value-of select="@Definite"/></td>
                                    </tr>
                                </table></div>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="Chromosome != '' and position() != 1">
                        <xsl:for-each select="Chromosome">
                            <xsl:if test="@Class='Single'">
                                <xsl:variable name="totalCost" select="translate(substring-before(@Cost,'-'),' ','')"/>
                                <xsl:variable name="rearrangementCost" select="@Rearrangment_Cost"/>
                                <h2><span>transformation map</span></h2>
                                <div class="transformations"><table>
                                    <tr>
                                        <th>character</th>
                                        <td><xsl:value-of select="translate(@Name,'_',' ')"/></td>
                                    </tr><tr>
                                        <th>type</th>
                                        <td>chromosome</td>
                                    </tr><tr>
                                        <th>edit cost</th>
                                        <td><xsl:value-of select="$totalCost - $rearrangementCost"/></td>
                                    </tr><tr>
                                        <th>rearrangement cost</th>
                                        <td><xsl:value-of select="$rearrangementCost"/></td>
                                    </tr><tr>
                                        <th>total cost</th>
                                        <td><xsl:value-of select="$totalCost"/></td>
                                    </tr><tr>
                                        <th>definite</th>
                                        <td><xsl:value-of select="@Definite"/></td>
                                    </tr><tr>
                                        <xsl:variable name="ancestorSequence" select="normalize-space(Sequence)"/>
                                        <th>ancestor sequence</th>
                                        <td><xsl:value-of select="$ancestorSequence"/></td>
                                    </tr><tr>
                                        <xsl:variable name="descendantReferenceCode" select="ChromosomeMap/SegmentMap/@DescendantReferenceCode"/>
                                        <xsl:variable name="descendantSequence" select="normalize-space(../..//Chromosome[@ReferenceCode = $descendantReferenceCode][@Class = 'Single']/Sequence)"/>
                                        <th>descendant sequence</th>
                                        <td><xsl:value-of select="$descendantSequence"/></td>
                                    </tr>
                                </table></div>
                                <table class="map">
                                    <tr>
                                        <th>from</th>
                                        <th>to</th>
                                        <th>type</th>
                                    </tr>
                                    <xsl:choose>
                                        <xsl:when test="ChromosomeMap/SegmentMap/@AncestorStartPosition != ''">
                                        <!-- non-annotated -->
                                            <xsl:for-each select="ChromosomeMap//SegmentMap">
                                                <!-- output ancestor set -->
                                                <tr>
                                                    <td>
                                                        <xsl:text>[</xsl:text>
                                                        <xsl:choose>
                                                            <xsl:when test="@AncestorDirection != '+'">
                                                                <xsl:value-of select="@AncestorEndPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@AncestorStartPosition"/>
                                                            </xsl:when>
                                                            <xsl:when test="@AncestorStartPosition = -1">
                                                                <xsl:text/>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:value-of select="@AncestorStartPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@AncestorEndPosition"/>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                        <xsl:text>]</xsl:text>
                                                    </td>
                                                    <!-- output descendant set -->
                                                    <td>
                                                        <xsl:text>[</xsl:text>
                                                        <xsl:choose>
                                                            <xsl:when test="@DescendantDirection != '+'">
                                                                <xsl:value-of select="@DescendantEndPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@DescendantStartPosition"/>
                                                            </xsl:when>
                                                            <xsl:when test="@DescendantStartPosition = -1">
                                                                <xsl:text/>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:value-of select="@DescendantStartPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@DescendantEndPosition"/>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                        <xsl:text>]</xsl:text>
                                                    </td>
                                                <!-- output type of trans -->
                                                    <td>
                                                        <xsl:choose>
                                                            <xsl:when test="@AncestorStartPosition = -1">
                                                                <xsl:text>origin</xsl:text>
                                                            </xsl:when>
                                                            <xsl:when test="@DescendantStartPosition = -1">
                                                                <xsl:text>loss</xsl:text>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:variable name="change" select="@AncestorStartPosition != @DescendantStartPosition or @AncestorEndPosition != @DescendantEndPosition"/>
                                                                <xsl:variable name="inversion" select="@AncestorDirection != @DescendantDirection"/>
                                                                <xsl:variable name="transformation" select="$change or $inversion"/>
                                                                <xsl:if test="$change">
                                                                    <xsl:text>change </xsl:text>
                                                                </xsl:if>
                                                                <xsl:if test="$inversion">
                                                                    <xsl:text>inversion</xsl:text>
                                                                </xsl:if>
                                                                <xsl:if test="not($transformation)">
                                                                    <xsl:text>no change</xsl:text>
                                                                </xsl:if>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                    </td>
                                                </tr>
                                            </xsl:for-each>
                                        </xsl:when>
                                        <xsl:otherwise>
                                        <!-- annotated -->
                                            <xsl:for-each select="ChromosomeMap//SegmentMap">
                                                <tr>
                                                    <td>
                                                        <xsl:text>[</xsl:text>
                                                        <xsl:value-of select="@AncestorSequenceOrder"/>
                                                        <xsl:text>]</xsl:text>
                                                    </td>
                                                    <td>
                                                        <xsl:text>[</xsl:text>
                                                        <xsl:value-of select="@DescendantSequenceOrder"/>
                                                        <xsl:text>]</xsl:text>
                                                    </td>
                                                    <td>
                                                        <xsl:choose>
                                                            <xsl:when test="@AncestorSequenceOrder = -1">
                                                                <xsl:text>origin</xsl:text>
                                                            </xsl:when>
                                                            <xsl:when test="@DescendantSequenceOrder = -1">
                                                                <xsl:text>loss</xsl:text>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:choose>
                                                                    <xsl:when test="@AncestorSequenceOrder != @DescendantSequenceOrder or @AncestorDirection != @DescendantDirection">
                                                                        <xsl:if test="@AncestorSequenceOrder != @DescendantSequenceOrder">
                                                                            <xsl:text>rearrangement </xsl:text>
                                                                        </xsl:if>
                                                                        <xsl:if test="@AncestorDirection != @DescendantDirection">
                                                                            <xsl:text>inversion</xsl:text>
                                                                        </xsl:if>
                                                                    </xsl:when>
                                                                    <xsl:otherwise>
                                                                        <xsl:text>no change</xsl:text>
                                                                    </xsl:otherwise>
                                                                </xsl:choose>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                    </td>
                                                </tr>
                                            </xsl:for-each>
                                        </xsl:otherwise>
                                    </xsl:choose>
                                </table>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="Genome != '' and position() != 1">
                        <xsl:for-each select="Genome">
                            <xsl:if test="@Class='Single'">
                                <xsl:variable name="totalCost" select="translate(substring-before(@Cost,'-'),' ','')"/>
                                <xsl:variable name="rearrangementCost" select="@Rearrangment_Cost"/>
                                <h2><span>transformation map</span></h2>
                                <div class="transformations"><table>
                                    <tr>
                                        <th>character</th>
                                        <td><xsl:value-of select="translate(@Name,'_',' ')"/></td>
                                    </tr><tr>
                                        <th>type</th>
                                        <td>multi-chromosome</td>
                                    </tr><tr>
                                        <th>edit cost</th>
                                        <td><xsl:value-of select="$totalCost - $rearrangementCost"/></td>
                                    </tr><tr>
                                        <th>rearrangement cost</th>
                                        <td><xsl:value-of select="$rearrangementCost"/></td>
                                    </tr><tr>
                                        <th>total cost</th>
                                        <td><xsl:value-of select="$totalCost"/></td>
                                    </tr><tr>
                                        <th>definite</th>
                                        <td><xsl:value-of select="@Definite"/></td>
                                    </tr><tr>
                                        <xsl:variable name="ancestorReferenceCode" select="GenomeMap/@AncestorReferenceCode"/>
                                        <xsl:variable name="ancestorSequence" select="normalize-space(../../../Node/Genome[@ReferenceCode = $ancestorReferenceCode][@Class='Single']/Sequence)"/>
                                        <th>ancestor sequence</th>
                                        <td><xsl:value-of select="$ancestorSequence"/></td>
                                    </tr><tr>
                                        <xsl:variable name="descendantSequence" select="normalize-space(Sequence)"/>
                                        <th>descendant sequence</th>
                                        <td><xsl:value-of select="$descendantSequence"/></td>
                                    </tr>
                                </table></div>
                                <table class="map">
                                    <tr>
                                        <th>from</th>
                                        <th>to</th>
                                        <th>type</th>
                                    </tr>
                                    <xsl:choose>
                                        <xsl:when test="GenomeMap/ChromosomeMap/SegmentMap/@AncestorStartPosition != ''">
                                        <!-- non-annotated -->
                                            <xsl:for-each select="GenomeMap/ChromosomeMap//SegmentMap">
                                                <tr>
                                                    <!-- output ancestor set -->
                                                    <td>
                                                        <xsl:text>c</xsl:text>
                                                        <xsl:value-of select="@AncestorChromosomeID"/>
                                                        <xsl:text>.:[</xsl:text>
                                                        <xsl:choose>
                                                            <xsl:when test="@AncestorDirection != '+'">
                                                                <xsl:value-of select="@AncestorEndPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@AncestorStartPosition"/>
                                                            </xsl:when>
                                                            <xsl:when test="@AncestorStartPosition = -1">
                                                                <xsl:text/>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:value-of select="@AncestorStartPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@AncestorEndPosition"/>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                        <xsl:text>]</xsl:text>
                                                    </td>
                                                    <!-- output descendant set -->
                                                    <td>
                                                        <xsl:text>c</xsl:text>
                                                        <xsl:value-of select="@DescendantChromosomeID"/>
                                                        <xsl:text>.:[</xsl:text>
                                                        <xsl:choose>
                                                            <xsl:when test="@DescendantDirection != '+'">
                                                                <xsl:value-of select="@DescendantEndPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@DescendantStartPosition"/>
                                                            </xsl:when>
                                                            <xsl:when test="@DescendantStartPosition = -1">
                                                                <xsl:text/>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:value-of select="@DescendantStartPosition"/>
                                                                <xsl:text>,</xsl:text>
                                                                <xsl:value-of select="@DescendantEndPosition"/>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                        <xsl:text>]</xsl:text>
                                                    </td>
                                                    <!-- output type of trans -->
                                                    <td>
                                                        <xsl:choose>
                                                            <xsl:when test="@AncestorStartPosition = -1">
                                                                <xsl:text>origin</xsl:text>
                                                            </xsl:when>
                                                            <xsl:when test="@DescendantStartPosition = -1">
                                                                <xsl:text>loss</xsl:text>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:variable name="change" select="@AncestorStartPosition != @DescendantStartPosition or @AncestorEndPosition != @DescendantEndPosition"/>
                                                                <xsl:variable name="inversion" select="@AncestorDirection != @DescendantDirection"/>
                                                                <xsl:variable name="jump" select="@AncestorChromosomeID != @DescendantChromosomeID"/>
                                                                <xsl:variable name="transformation" select="$change or $inversion or $jump"/>
                                                                <xsl:if test="$change">
                                                                    <xsl:text>change </xsl:text>
                                                                </xsl:if>
                                                                <xsl:if test="$inversion">
                                                                    <xsl:text>inversion </xsl:text>
                                                                </xsl:if>
                                                                <xsl:if test="$jump">
                                                                    <xsl:text>jump</xsl:text>
                                                                </xsl:if>
                                                                <xsl:if test="not($transformation)">
                                                                    <xsl:text>no change</xsl:text>
                                                                </xsl:if>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                    </td>
                                                </tr>
                                            </xsl:for-each>
                                        </xsl:when>
                                        <xsl:otherwise>
                                        <!-- annotated -->
                                            <xsl:for-each select="ChromosomeMap//SegmentMap">
                                                <tr>
                                                    <td>
                                                        <xsl:text>&#10;       [</xsl:text>
                                                        <xsl:value-of select="@AncestorSequenceOrder"/>
                                                        <xsl:text>]</xsl:text>
                                                    </td><td>
                                                        <xsl:text>[</xsl:text>
                                                        <xsl:value-of select="@DescendantSequenceOrder"/>
                                                        <xsl:text>]</xsl:text>
                                                    </td><td>
                                                        <xsl:choose>
                                                            <xsl:when test="@AncestorSequenceOrder = -1">
                                                                <xsl:text>origin</xsl:text>
                                                            </xsl:when>
                                                            <xsl:when test="@DescendantSequenceOrder = -1">
                                                                <xsl:text>loss</xsl:text>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <xsl:choose>
                                                                    <xsl:when test="@AncestorSequenceOrder != @DescendantSequenceOrder or @AncestorDirection != @DescendantDirection or @AncestorChromosomeID != @DescendantChromosomeID">
                                                                        <xsl:if test="@AncestorSequenceOrder != @DescendantSequenceOrder">
                                                                            <xsl:text>rearrangement </xsl:text>
                                                                        </xsl:if>
                                                                        <xsl:if test="@AncestorDirection != @DescendantDirection">
                                                                            <xsl:text>inversion </xsl:text>
                                                                        </xsl:if>
                                                                        <xsl:if test="@AncestorChromosomeID != @DescendantChromosomeID">
                                                                            <xsl:text>jump</xsl:text>
                                                                        </xsl:if>
                                                                    </xsl:when>
                                                                    <xsl:otherwise>
                                                                        <xsl:text>no change</xsl:text>
                                                                    </xsl:otherwise>
                                                                </xsl:choose>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                    </td>
                                                </tr>
                                            </xsl:for-each>
                                        </xsl:otherwise>
                                    </xsl:choose>
                                </table>
                            </xsl:if>
                        </xsl:for-each>
                    </xsl:if>
                    <xsl:if test="@Name != 'root'">
                        <p>ancestor : <xsl:value-of select="../../Node/@Name"/></p>
                    </xsl:if>
                </div>
            </xsl:for-each>
        </div>
    </xsl:template>
    
</xsl:stylesheet>
