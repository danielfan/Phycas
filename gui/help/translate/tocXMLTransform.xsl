<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0" xmlns:phyc="http://commandLanguage.phycas.org">
<xsl:output method="xml" indent="yes"/>
<xsl:strip-space elements="phyc:command_list"/>
<xsl:template match="/">

<toc version="1.0">
  <tocitem text="Phycas">
    <tocitem text="Introduction" target="intro"/>
    <tocitem text="Using the Phycas GUI" target="usingGui"/>
    <tocitem text="Commands" target="commands">
            <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:choose>
                    <xsl:when test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                        <xsl:for-each select="./command">
                            <xsl:sort select="@gui_label" data-type="number"/>
                            <xsl:call-template name="command_toc"/>
                        </xsl:for-each>
                    </xsl:when>
				</xsl:choose>
			</xsl:for-each>
    </tocitem>
  </tocitem>
</toc>

</xsl:template>

<xsl:template name="command_toc">
	<xsl:variable name="command_label"><xsl:value-of select="./@label"/></xsl:variable>
    <tocitem text="{$command_label}" target="command{$command_label}">
		<xsl:for-each select="cmd_param_group">
			<xsl:for-each select="cmd_param">
				<tocitem text="{./gui_description}" target="command{$command_label}Param{./@label}{./@placement}"/>
			</xsl:for-each>
		</xsl:for-each>
	</tocitem>
</xsl:template>


</xsl:stylesheet>
