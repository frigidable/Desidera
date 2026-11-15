import QtQuick 2.15
import QtQuick.Controls 2.15

// Isolated Observatory Studio comparison, no simulated science results.
// All project names, statuses and commands originate in the live Widgets hub.
Rectangle {
    id: root
    objectName: "observatory_qml_hub_root"
    width: 1180
    height: 740
    property bool darkMode: hubBridge.darkTheme
    property color backgroundColor: darkMode ? "#0C1424" : "#F2F5FA"
    property color surfaceColor: darkMode ? "#172238" : "#FFFFFF"
    property color raisedColor: darkMode ? "#1D2C46" : "#E9F0FA"
    property color textColor: darkMode ? "#E9F1FD" : "#1A2B45"
    property color mutedColor: darkMode ? "#AABAD1" : "#4D617A"
    property color borderColor: darkMode ? "#344862" : "#CED9E8"
    property color accentColor: darkMode ? "#65B8FA" : "#215CAA"
    property color selectionColor: darkMode ? "#1A3552" : "#E3F0FA"
    property color warnColor: darkMode ? "#EEC47A" : "#865313"
    property color warnBackground: darkMode ? "#33291D" : "#FFF3DE"
    color: backgroundColor

    Rectangle {
        id: topbar
        x: 22; y: 18; width: root.width - 44; height: 64
        radius: 12; color: root.surfaceColor
        border.width: 1; border.color: root.borderColor
        // Vector orbit mark instead of a font-dependent Unicode glyph.
        Rectangle {
            x: 20; y: 18; width: 27; height: 27; radius: 14
            color: "transparent"
            border.width: 2
            border.color: root.accentColor
            Rectangle {
                anchors.centerIn: parent
                width: 8; height: 8; radius: 4
                color: root.accentColor
            }
        }
        Text {
            x: 60; y: 15
            text: "DESIDERA"
            font.pixelSize: 21
            font.weight: Font.DemiBold
            color: root.textColor
        }
        Text {
            x: 224; y: 22
            text: "OBSERVATORY STUDIO    /    PROJECT HUB"
            font.pixelSize: 12
            font.letterSpacing: 0.8
            color: root.mutedColor
        }
        Rectangle {
            width: 154; height: 27; radius: 7
            anchors.right: parent.right
            anchors.rightMargin: 20
            anchors.verticalCenter: parent.verticalCenter
            color: root.raisedColor
            Text {
                anchors.centerIn: parent
                text: "QML comparison · beta"
                font.pixelSize: 11
                color: root.mutedColor
            }
        }
    }

    Rectangle {
        id: library
        x: 22; y: 98
        width: Math.max(338, Math.min(440, root.width * 0.36))
        height: root.height - y - 56
        radius: 12
        color: root.surfaceColor
        border.width: 1
        border.color: root.borderColor

        Text {
            id: libraryHeading
            x: 20; y: 19
            text: "Galaxy projects"
            font.pixelSize: 21
            font.weight: Font.DemiBold
            color: root.textColor
        }
        Text {
            x: 20; y: 54
            text: hubBridge.visibleCount === hubBridge.projectCount
                ? hubBridge.projectCount + " saved galaxy projects"
                : hubBridge.visibleCount + " of " + hubBridge.projectCount + " galaxies"
            font.pixelSize: 12
            color: root.mutedColor
        }
        TextField {
            id: queryField
            objectName: "qml_hub_filter"
            x: 20; y: 87; width: library.width - 40; height: 40
            placeholderText: "Search galaxies or paths..."
            font.pixelSize: 13
            color: root.textColor
            placeholderTextColor: root.mutedColor
            selectByMouse: true
            onTextChanged: hubBridge.setFilter(text)
            background: Rectangle {
                radius: 8
                color: root.raisedColor
                border.width: queryField.activeFocus ? 2 : 1
                border.color: queryField.activeFocus ? root.accentColor : root.borderColor
            }
        }
        Text {
            x: 20; y: 144
            text: "RECENT GALAXIES"
            color: root.accentColor
            font.pixelSize: 11
            font.weight: Font.DemiBold
            font.letterSpacing: 1
        }
        ListView {
            id: galaxyList
            objectName: "qml_hub_list"
            x: 12; y: 169
            width: library.width - 24
            height: library.height - 270
            clip: true
            spacing: 8
            focus: true
            model: projectModel
            // Never commit a different project on model reset or theme edit.
            keyNavigationEnabled: false
            Component.onCompleted: currentIndex = hubBridge.selectedIndex
            Connections {
                target: hubBridge
                function onChanged() {
                    galaxyList.currentIndex = hubBridge.selectedIndex
                }
            }
            Keys.onDownPressed: {
                if (galaxyList.count > 0) {
                    galaxyList.currentIndex = Math.min(
                        galaxyList.count - 1, galaxyList.currentIndex + 1)
                    hubBridge.choose(hubBridge.idAt(galaxyList.currentIndex))
                }
                event.accepted = true
            }
            Keys.onUpPressed: {
                if (galaxyList.count > 0) {
                    galaxyList.currentIndex = Math.max(
                        0, galaxyList.currentIndex - 1)
                    hubBridge.choose(hubBridge.idAt(galaxyList.currentIndex))
                }
                event.accepted = true
            }
            Keys.onReturnPressed: hubBridge.openSelected()
            Keys.onEnterPressed: hubBridge.openSelected()
            delegate: Rectangle {
                required property int index
                required property string galaxyName
                required property string statusText
                required property string coordinates
                required property string projectId
                required property bool isAvailable
                width: galaxyList.width
                height: 93
                radius: 9
                color: hubBridge.selectedId === projectId
                    ? root.selectionColor : root.raisedColor
                border.width: hubBridge.selectedId === projectId ? 2 : 1
                border.color: hubBridge.selectedId === projectId
                    ? root.accentColor : root.borderColor
                Rectangle {
                    x: 13; y: 16; width: 40; height: 40; radius: 9
                    color: root.darkMode ? "#253955" : "#D9E8F6"
                    Rectangle {
                        anchors.centerIn: parent
                        width: 24; height: 24; radius: 12
                        color: "transparent"
                        border.width: 2
                        border.color: root.accentColor
                        Rectangle {
                            anchors.centerIn: parent
                            width: 7; height: 7; radius: 4
                            color: root.accentColor
                        }
                    }
                }
                Text {
                    x: 64; y: 13
                    width: parent.width - 79
                    text: galaxyName
                    elide: Text.ElideRight
                    color: root.textColor
                    font.pixelSize: 15
                    font.weight: Font.DemiBold
                }
                Text {
                    x: 64; y: 39
                    width: parent.width - 79
                    text: statusText
                    elide: Text.ElideRight
                    color: isAvailable ? root.accentColor : root.warnColor
                    font.pixelSize: 12
                }
                Text {
                    x: 64; y: 63
                    width: parent.width - 79
                    text: coordinates
                    elide: Text.ElideRight
                    color: root.mutedColor
                    font.pixelSize: 11
                }
                MouseArea {
                    anchors.fill: parent
                    onClicked: hubBridge.choose(projectId)
                    onDoubleClicked: {
                        hubBridge.choose(projectId)
                        hubBridge.openSelected()
                    }
                }
            }
            ScrollBar.vertical: ScrollBar { policy: ScrollBar.AsNeeded }
        }
        Text {
            anchors.centerIn: galaxyList
            width: galaxyList.width - 44
            visible: hubBridge.visibleCount === 0
            text: hubBridge.projectCount === 0
                ? "No saved galaxies yet.\nCreate a new project or import a folder."
                : "No matching galaxies.\nTry another search."
            wrapMode: Text.WordWrap
            horizontalAlignment: Text.AlignHCenter
            lineHeight: 1.4
            color: root.mutedColor
            font.pixelSize: 14
        }
        Button {
            id: newProjectButton
            objectName: "qml_hub_new_project"
            x: 20; y: library.height - 85
            width: (library.width - 48) * 0.53
            height: 45
            text: "+  New galaxy"
            onClicked: hubBridge.createProject()
            background: Rectangle {
                radius: 8
                color: newProjectButton.down ? "#16538E" : root.accentColor
            }
            contentItem: Text {
                text: newProjectButton.text
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
                font.pixelSize: 13
                font.weight: Font.DemiBold
                color: root.darkMode ? "#071629" : "#FFFFFF"
            }
        }
        Button {
            id: importButton
            objectName: "qml_hub_import_folder"
            x: newProjectButton.x + newProjectButton.width + 8
            y: newProjectButton.y
            width: library.width - x - 20
            height: 45
            text: "Import folder"
            onClicked: hubBridge.importFolder()
            background: Rectangle {
                radius: 8
                color: root.raisedColor
                border.width: 1
                border.color: root.borderColor
            }
            contentItem: Text {
                text: importButton.text
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
                font.pixelSize: 13
                color: root.textColor
            }
        }
    }

    Rectangle {
        id: inspector
        x: library.x + library.width + 16
        y: library.y
        width: root.width - x - 22
        height: library.height
        radius: 12
        color: root.surfaceColor
        border.width: 1
        border.color: root.borderColor

        Rectangle {
            x: 20; y: 20; width: 76; height: 76; radius: 11
            color: root.raisedColor
            border.width: 1
            border.color: root.borderColor
            Image {
                id: cachedGalaxyPreview
                anchors.fill: parent
                anchors.margins: 2
                source: hubBridge.previewUrl
                sourceSize.width: 180
                sourceSize.height: 180
                fillMode: Image.PreserveAspectCrop
                asynchronous: true
                visible: status === Image.Ready
            }
            Rectangle {
                visible: !cachedGalaxyPreview.visible
                anchors.centerIn: parent
                width: 44; height: 44; radius: 22
                color: "transparent"
                border.width: 3
                border.color: root.accentColor
                Rectangle {
                    anchors.centerIn: parent
                    width: 12; height: 12; radius: 6
                    color: root.accentColor
                }
            }
        }
        Text {
            x: 114; y: 20
            width: inspector.width - 260
            text: hubBridge.selectedName
            elide: Text.ElideRight
            font.pixelSize: 25
            font.weight: Font.DemiBold
            color: root.textColor
        }
        Text {
            x: 114; y: 61
            width: inspector.width - 146
            text: hubBridge.selectedCoords
            wrapMode: Text.WordWrap
            font.pixelSize: 12
            color: root.mutedColor
        }
        Rectangle {
            x: inspector.width - 144; y: 21
            width: 123; height: 27; radius: 7
            color: root.warnBackground
            Text {
                anchors.centerIn: parent
                text: "SCIENCE HOLD"
                color: root.warnColor
                font.pixelSize: 12
                font.weight: Font.DemiBold
            }
        }
        Text {
            x: 114; y: 100
            width: inspector.width - 134
            text: hubBridge.previewLabel
            textFormat: Text.PlainText
            color: root.mutedColor
            font.pixelSize: 10
            elide: Text.ElideRight
        }

        property bool detailsExpanded: false

        // Scroll evidence, keep the action toolbar and primary Open pinned.
        ScrollView {
            id: evidenceScroll
            objectName: "qml_hub_evidence_scroll"
            x: 20; y: 119
            width: inspector.width - 40
            height: Math.max(100, inspector.height - 346)
            clip: true
            ScrollBar.vertical: ScrollBar {
                policy: ScrollBar.AsNeeded
                width: 9
                background: Rectangle { color: "transparent" }
                contentItem: Rectangle {
                    implicitWidth: 7
                    radius: 3
                    color: root.darkMode ? "#466282" : "#AFBFCE"
                }
            }
            Column {
                id: evidenceColumn
                width: evidenceScroll.availableWidth - 2
                spacing: 6

                Text {
                    text: "SAVED PROJECT EVIDENCE"
                    color: root.accentColor
                    font.pixelSize: 11
                    font.weight: Font.DemiBold
                    font.letterSpacing: 0.9
                }
                Text {
                    width: parent.width
                    text: "Preparation:  " + hubBridge.selectedStatus
                    textFormat: Text.PlainText
                    color: root.textColor
                    font.pixelSize: 13
                    font.weight: Font.DemiBold
                    wrapMode: Text.WordWrap
                }
                Text {
                    width: parent.width
                    text: "Dataset:  " + hubBridge.datasetLabel
                    textFormat: Text.PlainText
                    color: root.mutedColor
                    font.pixelSize: 12
                    wrapMode: Text.WordWrap
                }
                Text {
                    width: parent.width
                    text: "FITS:  " + hubBridge.fitsLabel
                    textFormat: Text.PlainText
                    color: root.mutedColor
                    font.pixelSize: 12
                    wrapMode: Text.WordWrap
                }
                Rectangle {
                    width: parent.width; height: 1
                    color: root.borderColor
                }
                Text {
                    text: "CONFIRMED CATALOGUE IDENTITIES"
                    font.pixelSize: 11
                    font.letterSpacing: 0.6
                    font.weight: Font.DemiBold
                    color: root.accentColor
                }
                Text {
                    objectName: "qml_catalogue_summary"
                    width: parent.width
                    text: hubBridge.catalogueLabel
                    textFormat: Text.PlainText
                    wrapMode: Text.WordWrap
                    color: root.textColor
                    font.pixelSize: 13
                }
                Text {
                    width: parent.width
                    text: hubBridge.literatureLabel
                    visible: text.length > 0
                    textFormat: Text.PlainText
                    wrapMode: Text.WordWrap
                    color: root.mutedColor
                    font.pixelSize: 12
                }
                Rectangle {
                    width: parent.width; height: 1
                    color: root.borderColor
                }
                Text {
                    text: "FITS / WCS QA"
                    font.pixelSize: 11
                    font.letterSpacing: 0.6
                    font.weight: Font.DemiBold
                    color: root.accentColor
                }
                Text {
                    objectName: "qml_qa_summary"
                    width: parent.width
                    text: hubBridge.qaSummary
                    textFormat: Text.PlainText
                    wrapMode: Text.WordWrap
                    color: root.textColor
                    font.pixelSize: 13
                }
                Text {
                    width: parent.width
                    text: "Previous diagnostics are evidence, not accepted measurements."
                    color: root.mutedColor
                    font.pixelSize: 11
                    wrapMode: Text.WordWrap
                }
                Button {
                    id: evidenceToggle
                    objectName: "qml_hub_evidence_toggle"
                    width: parent.width
                    height: 34
                    enabled: hubBridge.selectedId.length > 0
                    text: inspector.detailsExpanded
                        ? "Hide sources and QA details"
                        : "Show sources and QA details"
                    onClicked: inspector.detailsExpanded = !inspector.detailsExpanded
                    background: Rectangle {
                        radius: 8
                        color: root.raisedColor
                        border.width: 1
                        border.color: root.borderColor
                    }
                    contentItem: Text {
                        text: evidenceToggle.text
                        font.pixelSize: 12
                        color: root.accentColor
                        verticalAlignment: Text.AlignVCenter
                        horizontalAlignment: Text.AlignHCenter
                    }
                }
                Text {
                    objectName: "qml_hub_catalogue_sources"
                    visible: inspector.detailsExpanded && hubBridge.selectedId.length > 0
                    width: parent.width
                    text: "CATALOGUE SOURCES\n" + hubBridge.catalogueEvidence
                    textFormat: Text.PlainText
                    wrapMode: Text.WordWrap
                    color: root.mutedColor
                    font.pixelSize: 12
                }
                Text {
                    objectName: "qml_hub_qa_sources"
                    visible: inspector.detailsExpanded && hubBridge.selectedId.length > 0
                    width: parent.width
                    text: "FITS/WCS AUDIT EVIDENCE\n" + hubBridge.qaEvidence
                    textFormat: Text.PlainText
                    wrapMode: Text.WordWrap
                    color: root.mutedColor
                    font.pixelSize: 12
                }
                Text {
                    objectName: "qml_hub_project_sources"
                    visible: inspector.detailsExpanded && hubBridge.selectedId.length > 0
                    width: parent.width
                    text: "PROJECT INFORMATION\n" + hubBridge.projectEvidence
                    textFormat: Text.PlainText
                    wrapMode: Text.WordWrap
                    color: root.mutedColor
                    font.pixelSize: 12
                }
            }
        }

        Rectangle {
            id: nextBox
            objectName: "qml_hub_next_action_card"
            x: 20; y: inspector.height - 216
            width: inspector.width - 40; height: 100
            radius: 9
            color: root.raisedColor
            border.width: 1; border.color: root.borderColor
            Text {
                x: 13; y: 10
                text: "RECOMMENDED NEXT STEP"
                color: root.accentColor
                font.pixelSize: 11
                font.weight: Font.DemiBold
            }
            Text {
                x: 13; y: 32
                width: nextBox.width - 26
                height: 44
                text: hubBridge.nextAction
                color: root.textColor
                textFormat: Text.PlainText
                font.pixelSize: 13
                font.weight: Font.DemiBold
                wrapMode: Text.WordWrap
                elide: Text.ElideRight
            }
            Text {
                x: 13; y: 80
                width: nextBox.width - 26
                text: "From project readiness · no science adopted"
                color: root.mutedColor
                font.pixelSize: 10
                elide: Text.ElideRight
            }
        }

        Button {
            id: catalogueButton
            objectName: "qml_hub_catalogue_actions"
            x: 20; y: inspector.height - 106
            width: (inspector.width - 49) / 2; height: 35
            text: "Catalogues"
            enabled: hubBridge.canCatalogue
            onClicked: catalogueMenu.open()
            background: Rectangle {
                radius: 8; color: root.raisedColor
                border.width: 1; border.color: root.borderColor
            }
            contentItem: Text {
                text: catalogueButton.text
                color: catalogueButton.enabled ? root.textColor : root.mutedColor
                font.pixelSize: 12
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
            }
            Menu {
                id: catalogueMenu
                y: parent.height
                MenuItem {
                    text: "SIMBAD · inspect identity…"
                    onTriggered: hubBridge.runTool("simbad")
                }
                MenuItem {
                    text: "NED · inspect identity…"
                    onTriggered: hubBridge.runTool("ned")
                }
            }
        }
        Button {
            id: projectToolsButton
            objectName: "qml_hub_project_tools"
            x: catalogueButton.x + catalogueButton.width + 9
            y: catalogueButton.y
            width: inspector.width - x - 20; height: 35
            text: "Project tools"
            enabled: hubBridge.canRemove
            onClicked: projectToolsMenu.open()
            background: Rectangle {
                radius: 8; color: root.raisedColor
                border.width: 1; border.color: root.borderColor
            }
            contentItem: Text {
                text: projectToolsButton.text
                color: projectToolsButton.enabled ? root.textColor : root.mutedColor
                font.pixelSize: 12
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
            }
            Menu {
                id: projectToolsMenu
                y: parent.height
                MenuItem {
                    text: "Analysis configurations…"
                    enabled: hubBridge.canConfigure
                    onTriggered: hubBridge.runTool("config")
                }
                MenuItem {
                    text: "Regenerate cached FITS thumbnail"
                    enabled: hubBridge.canLocalPreview
                    onTriggered: hubBridge.runTool("localPreview")
                }
                MenuItem {
                    text: "Fetch colour preview (network request)"
                    enabled: hubBridge.canRemotePreview
                    onTriggered: hubBridge.runTool("remotePreview")
                }
                MenuSeparator {}
                MenuItem {
                    text: "Remove from recent list…"
                    onTriggered: hubBridge.runTool("remove")
                }
            }
        }
        Button {
            id: openButton
            objectName: "qml_hub_open_project"
            x: 20; y: inspector.height - 62
            width: inspector.width - 40; height: 44
            text: hubBridge.openLabel
            enabled: hubBridge.selectedId.length !== 0
            onClicked: hubBridge.openSelected()
            background: Rectangle {
                radius: 8
                color: openButton.enabled ? root.accentColor : root.raisedColor
            }
            contentItem: Text {
                text: openButton.text
                color: openButton.enabled
                    ? (root.darkMode ? "#071629" : "#FFFFFF") : root.mutedColor
                font.pixelSize: 14
                font.weight: Font.DemiBold
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
            }
        }
    }

    Text {
        x: 28; y: root.height - 34
        text: "Original FITS are read-only · catalogue associations are not scientific adoption"
        color: root.mutedColor
        font.pixelSize: 11
    }
    Text {
        anchors.right: parent.right
        anchors.rightMargin: 28
        y: root.height - 34
        text: "Qt Quick experimental preview"
        color: root.mutedColor
        font.pixelSize: 11
    }
}
