<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>
<xsl:strip-space elements="*"/>

<xsl:key name="alphabet" match="Element" use="@Value" />

	<xsl:template match="/">
		<xsl:apply-templates select="Diagnosis/Data/Taxa"/>
		<xsl:apply-templates select="Diagnosis/Data/Characters"/>
	</xsl:template>

	<xsl:template match="Diagnosis/Data/Taxa">
		<xsl:text>BEGIN TAXA;</xsl:text>
		<xsl:text>&#10;&#x09;DIMENSIONS NTAX=</xsl:text>
		<xsl:value-of select="count(Taxon)"/>
		<xsl:text>;</xsl:text>
		<xsl:text>&#10;&#x09;TAXLABELS </xsl:text>
		<xsl:for-each select="Taxon">
			<xsl:value-of select="@Name"/>
			<xsl:text> </xsl:text>
		</xsl:for-each>
		<xsl:text>;</xsl:text>
		<xsl:text>&#10;END;</xsl:text>
	</xsl:template>
	
	<xsl:template match="Diagnosis/Data/Characters">
		<xsl:text>&#10;BEGIN CHARACTERS;</xsl:text>
		<xsl:text>&#10;&#x09;DIMENSIONS NCHAR=</xsl:text>
		<xsl:value-of select="count(Character)"/>
		<xsl:text>;</xsl:text>
		<xsl:text>&#10;&#x09;FORMAT DATATYPE=DNA MISSING=</xsl:text>
		<xsl:value-of select="Character/child::*/@Missing_Symbol"/>
		<xsl:text> GAP=</xsl:text>
		<xsl:value-of select="Character/child::*/@Gap_Symbol"/>
		
		<xsl:text> EQUATE="</xsl:text>
		<xsl:value-of select="Character/child::*/Equivalencies/Equivalency/@From"/>
		<xsl:text>=</xsl:text>
		<xsl:value-of select="Character/child::*/Equivalencies/Equivalency/@To"/>
		<xsl:text>"</xsl:text>
		
		<xsl:text> MATCHCHAR=</xsl:text>
		<xsl:value-of select="Character/child::*/@Match_Symbol"/>
		<xsl:text> ELIMINATE {</xsl:text>
		<xsl:for-each select="Character/child::*">
			<xsl:if test="@Ignore = 'true'">
				<xsl:value-of select="@Name"/>
				<xsl:text> </xsl:text>
			</xsl:if>
		</xsl:for-each>
		<xsl:text>} CHARSTATELABELS </xsl:text>
		<xsl:for-each select="Character/child::*">
			<xsl:value-of select="position()"/>
			<xsl:text> /{</xsl:text>
			<xsl:for-each select="Alphabet/Element">
				<xsl:value-of select="@Value"/>
				<xsl:text> </xsl:text>
			</xsl:for-each>
			<xsl:text>} </xsl:text>
		</xsl:for-each>
		<xsl:text>;</xsl:text>
		<xsl:text>&#10;END;</xsl:text>
	</xsl:template>
	
</xsl:stylesheet>
