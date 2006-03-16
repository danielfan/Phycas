<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0" xmlns:phyc="http://commandLanguage.phycas.org">
<xsl:output method="xml" indent="yes"/>
<xsl:strip-space elements="phyc:command_list"/>
<xsl:template match="/">

<index version="1.0">
    <indexitem text="Phycas" target="intro"/>
    <indexitem text="GUI" target="usingGui"/>
            <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:choose>
                    <xsl:when test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                        <xsl:for-each select="./command">
                            <xsl:sort select="@gui_label" data-type="number"/>
                            <xsl:call-template name="command_index"/>
                        </xsl:for-each>
                    </xsl:when>
				</xsl:choose>
			</xsl:for-each>
</index>

</xsl:template>

<xsl:template name="command_index">
	<xsl:variable name="command_label"><xsl:value-of select="./@label"/></xsl:variable>
    <indexitem text="{$command_label}" target="command{$command_label}"/>
</xsl:template>


</xsl:stylesheet>
