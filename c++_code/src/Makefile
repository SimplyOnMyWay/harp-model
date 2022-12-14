### libstk Makefile - for various flavors of unix

LIBRARY = libstk.a
AR = ar -rsc
RM = /bin/rm
OBJECT_PATH = Release
vpath %.o $(OBJECT_PATH)

OBJECTS	=	Stk.o Generator.o Noise.o SubNoise.o Blit.o BlitSaw.o BlitSquare.o \
					Envelope.o ADSR.o Asymp.o Modulate.o SingWave.o SineWave.o Granulate.o \
					FileRead.o FileWrite.o WvIn.o FileWvIn.o WaveLoop.o WvOut.o FileWvOut.o \
					Filter.o OneZero.o OnePole.o PoleZero.o TwoZero.o TwoPole.o \
					BiQuad.o FormSwep.o Delay.o DelayL.o DelayA.o \
					Effect.o PRCRev.o JCRev.o NRev.o \
					Chorus.o Echo.o PitShift.o \
					Function.o ReedTable.o JetTable.o BowTable.o \
					Voicer.o Vector3D.o Sphere.o \
					\
					Instrmnt.o Clarinet.o BlowHole.o Saxofony.o Flute.o Brass.o BlowBotl.o \
					Bowed.o Plucked.o StifKarp.o Sitar.o PluckTwo.o Mandolin.o Mesh2D.o \
					FM.o Rhodey.o Wurley.o TubeBell.o HevyMetl.o PercFlut.o BeeThree.o FMVoices.o \
					Sampler.o Moog.o Simple.o Drummer.o Shakers.o \
					Modal.o ModalBar.o BandedWG.o Resonate.o VoicForm.o Phonemes.o Whistle.o \
					\
					Messager.o Skini.o MidiFileIn.o


INCLUDE = 
ifeq ($(strip $(INCLUDE)),)
	vpath %.h ../include
	INCLUDE = -I../include
else
	vpath %.h ../include $(INCLUDE)
	INCLUDE = -I../include -I
endif

CC       = g++
DEFS     = -D__LITTLE_ENDIAN__
DEFS    += 
CFLAGS   = -O3
CFLAGS  += -Wall -g -D__GXX__ $(INCLUDE)

REALTIME = yes
ifeq ($(REALTIME),yes)
	OBJECTS += RtMidi.o RtAudio.o RtWvOut.o RtWvIn.o RtDuplex.o InetWvOut.o InetWvIn.o Thread.o Mutex.o Socket.o TcpClient.o TcpServer.o UdpSocket.o
	DEFS    += -D__LINUX_ALSA__ -D__LINUX_ALSASEQ__
endif

RAWWAVES = 
ifeq ($(strip $(RAWWAVES)), )
	RAWWAVES = ../rawwaves/
endif
DEFS    += -DRAWWAVE_PATH=\"$(RAWWAVES)\"

%.o : %.cpp
	$(CC) $(CFLAGS) $(DEFS) -c $(<) -o $(OBJECT_PATH)/$@

all : $(LIBRARY)

$(LIBRARY) : $(OBJECTS)
	/bin/rm -f $(LIBRARY)
	$(AR) $(LIBRARY) $(OBJECT_PATH)/*.o

$(OBJECTS) : Stk.h

clean : 
	-rm $(OBJECT_PATH)/*.o
	-rm $(LIBRARY)
