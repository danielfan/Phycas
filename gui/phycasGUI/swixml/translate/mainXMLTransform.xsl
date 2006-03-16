<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0" xmlns:phyc="http://commandLanguage.phycas.org">
<xsl:output method="xml" indent="yes"/>
<xsl:strip-space elements="phyc:command_list"/>
<xsl:template match="/">
<frame id="frame" size="800,600" layout="BorderLayout" title="P H Y C A S" DefaultCloseOperation="JFrame.EXIT_ON_CLOSE" plaf="com.sun.java.swing.plaf.windows.WindowsLookAndFeel">
  <menubar id="mb">
    <menu id="menu_main" text="FILE" mnemonic="VK_F">
      <menuitem id="mi_welcome" text="Welcome Screen" Action="showAction" ActionCommand="welcome" mnemonic="VK_W"/>
      <menuitem id="mi_exit" text="Exit" Action="exitAction" mnemonic="VK_E" macos_quit="true" />
    </menu>
            <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:choose>
                    <xsl:when test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
					<menu id="menu_{./@name}" text="{translate(./@name, 'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ')}" mnemonic="VK_{substring(translate(./@name, 'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'),1,1)}"> 
                        <xsl:for-each select="./command">
                            <xsl:sort select="@gui_order" data-type="number"/>
                            <xsl:call-template name="command_menu"/>
                        </xsl:for-each>
					</menu>
                    </xsl:when>
                    <xsl:otherwise><xsl:for-each select="./command">
                            <xsl:if test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                                <xsl:call-template name="command_menu"/>
                            </xsl:if>
                        </xsl:for-each>
                   </xsl:otherwise>
                </xsl:choose>
             </xsl:for-each>
    <menu id="menu_settings" text="SETTINGS" mnemonic="VK_S">
	<buttongroup>
		<checkboxmenuitem id="mi_settings_displayWarnings" text="Display Warnings" Action="action_settings_displayWarnings" ActionCommand="action_settings_displayWarnings" mnemonic="VK_D" selected="true"/>
		<checkboxmenuitem id="mi_settings_logWarnings" text="Log Warnings" Action="action_settings_logWarnings" ActionCommand="action_settings_logWarnings" mnemonic="VK_L"/>
	</buttongroup>
	<separator/>
    </menu>
    <menu id="menu_window" text="WINDOW" mnemonic="VK_W">
      <menuitem id="mi_output" text="Output" Action="showAction_Output" ActionCommand="showAction_Output" mnemonic="VK_O"/>
    </menu>
    <menu id="menu_help" text="HELP" mnemonic="VK_H">
      <menuitem id="mi_help" text="Launch Help" mnemonic="VK_H"/>
      <menuitem id="mi_help_context" text="Get Help on a Specific Item" mnemonic="VK_G"/>
    </menu>
  </menubar>

<panel id="pnl_main" Layout="BorderLayout" focusable="false">

<splitpane id="topPane" orientation="HORIZONTAL" resizeWeight="0.9" dividerLocation="0.25" focusable="false">

<scrollpane focusable="false" id="pnlScroll">
      <panel id="pnl" constraints="BorderLayout.CENTER" Layout="CardLayout(10,10)" Border="EtchedBorder" focusable="false">
  
    <panel id="pnl_welcome" Name="welcome" constraints="welcome" focusable="false">
        <vbox alignment="0.5" focusable="false">
			<hbox alignment="0.5" focusable="false">
				<label text="Welcome to Phycas!" Font="Serif-BOLD-18" Foreground="blue" focusable="false"/>
			</hbox>
			<hbox alignment="0.5" focusable="false">
				<label text="" focusable="false" maximumSize="20,20" preferredSize="20,20"/>
			</hbox>
			<hbox alignment="0.5" focusable="false">
				<label text="Please choose a command from one of the menus or enter a command below." Font="Serif-14" Foreground="black" focusable="false" maximumSize="450,20" preferredSize="450,20"/>
			</hbox>
			<hbox alignment="0.5" focusable="false">
				<label text="" focusable="false" maximumSize="10,10" preferredSize="10,10"/>
			</hbox>
			<hbox alignment="0.5" focusable="false">
				<textfield id="tfCommand" name="tfCommand" text="" alignmentY="0.0" maximumSize="450,20" preferredSize="450,20" toolTipText = "Manually enter a command."/>
			</hbox>
			<hbox alignment="0.5" focusable="false">
				<button text="OK" Action="okAction_Welcome" ActionCommand = "okAction_Welcome"/>
			</hbox>
       </vbox>
    </panel>
            <xsl:for-each select="//phyc:command_list/command_set">
                <xsl:if test="(./@user_interface='all') or (./@user_interface='gui') or (./@user_interface='batch|gui') or (./@user_interface='console|gui')">
                        <xsl:for-each select="./command">
                            <xsl:call-template name="command_panel"/></xsl:for-each>
                </xsl:if>

            </xsl:for-each>
  </panel>
</scrollpane>
<scrollpane focusable="false" id="taMsgsScroll">
        <textarea id="taMsgs" text="Welcome to Phycas!" rows="5" columns="50" autoscrolls="true" editable="false" linewrap="true"  wrapstyleword ="true" Font="Serif-14" focusable="false"/>
 </scrollpane>
</splitpane>
</panel>
</frame>
</xsl:template>

<xsl:template name="command_menu">
        <menuitem id="mi_{./@label}" text="{./@label}" action="showAction" actionCommand="pnl_{./@label}" mnemonic="VK_{substring(translate(./@label, 'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'),1,1)}"> </menuitem>
</xsl:template>

<xsl:template name="command_panel">
    <panel id="pnl_{./@label}" constraints="pnl_{./@label}" include="xml/Commands.xml#PNL_{./@label}"/>
</xsl:template>

</xsl:stylesheet>
