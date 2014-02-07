/*=========================================================================

   Program: ParaView
   Module:    pqPCLOpenNISource.h

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2. 

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "pqObjectPanel.h"

#include "pqView.h"

#include <vtkSMSourceProxy.h>
#include <vtkSMIntVectorProperty.h>
 
#include <QPushButton>
#include <QVBoxLayout>
#include <QCheckBox>
#include <QSlider>
#include <QLabel>
#include <QTimer>

class pqPCLOpenNISource : public pqObjectPanel
{

  Q_OBJECT

public:

  pqPCLOpenNISource(pqProxy* proxy, QWidget* p) : pqObjectPanel(proxy, p)
  {

    QVBoxLayout* layout = new QVBoxLayout(this);

    this->Timer = new QTimer(this);
    this->connect(this->Timer, SIGNAL(timeout()), SLOT(onPollSource()));

    QCheckBox* check = new QCheckBox("Enable OpenNI grabber thread");
    this->connect(check, SIGNAL(toggled(bool)), SLOT(onEnableOpenNIGrabberChecked(bool)));
    layout->addWidget(check);

    layout->addWidget(new QLabel);

    QPushButton* refreshButton = new QPushButton("Refresh");
    this->connect(refreshButton, SIGNAL(clicked()), SLOT(onRefreshClicked()));
    layout->addWidget(refreshButton);

    QCheckBox* autoRefreshCheck = new QCheckBox("Auto refresh");
    this->connect(autoRefreshCheck, SIGNAL(toggled(bool)), SLOT(onAutoRefreshChecked(bool)));
    layout->addWidget(autoRefreshCheck);

    QSlider* slider = new QSlider(Qt::Horizontal);
    this->connect(slider, SIGNAL(valueChanged(int)), SLOT(onSliderValueChanged(int)));
    layout->addWidget(slider);
    this->AutoRefreshLabel = new QLabel();
    layout->addWidget(this->AutoRefreshLabel);

    layout->addStretch();

    slider->setMinimum(1);
    slider->setMaximum(200);
    slider->setValue(100);
  }

public slots:

  void onEnableOpenNIGrabberChecked(bool checked)
  {
    vtkSMSourceProxy* sourceProxy = vtkSMSourceProxy::SafeDownCast(this->proxy());
    if (!sourceProxy)
      {
      return;
      }

    if (checked)
      {
      sourceProxy->InvokeCommand("StartGrabber");
      }
    else
      {
      sourceProxy->InvokeCommand("StopGrabber");
      }
  }

  void onAutoRefreshChecked(bool checked)
  {
    if (checked)
      {
      this->Timer->start();
      }
    else
      {
      this->Timer->stop();
      }
  }

  void onSliderValueChanged(int sliderValue)
  {
    int timeoutMax = 5000;
    int timeout = timeoutMax * sliderValue / 200.0;
    this->Timer->setInterval(timeout);
    this->AutoRefreshLabel->setText(QString("Auto refresh timeout: %0 s").arg(timeout/1000.0,  0, 'f', 2));
  }

  void onRefreshClicked()
  {
    this->onPollSource();
  }

  void onPollSource()
  {
    vtkSMSourceProxy* sourceProxy = vtkSMSourceProxy::SafeDownCast(this->proxy());
    if (!sourceProxy)
      {
      return;
      }

    vtkSMIntVectorProperty* prop = vtkSMIntVectorProperty::SafeDownCast(sourceProxy->GetProperty("HasNewData"));
    if (!prop)
      {
      return;
      }

    sourceProxy->InvokeCommand("Poll");
    sourceProxy->UpdatePropertyInformation(prop);
    int hasNewData = prop->GetElement(0);

    if (hasNewData && this->view())
      {
      this->view()->render();
      }
  }


protected:

  QLabel* AutoRefreshLabel;
  QTimer* Timer;
};
