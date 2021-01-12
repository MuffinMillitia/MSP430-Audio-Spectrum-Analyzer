#include <msp430.h> 
#include <string.h>
#include "QFFT.h"

/**
 * main.c
 */
#define ADC_BUFFER_SIZE 16                                  // Number of samples used for FFT
#define UART_BUFFER_SIZE 10                                 // Height of spectrum plot
#define SAMPLE_FREQUENCY 62500                              // MUST BE COMPUTED USING CLOCK FREQUENCY AND PRESCALAR ADC CONFIGURATION

static unsigned int ADC_buffer[ADC_BUFFER_SIZE];            // Stores ADC values
static unsigned volatile int ADC_index;
static volatile int ADC_buffer_flag;

static unsigned char UART_buffer[UART_BUFFER_SIZE][ADC_BUFFER_SIZE/2+1];    // Stores UART tx values
static unsigned volatile int UART_index;
static volatile int UART_buffer_flag = 0;

static unsigned volatile int spectrum_scale = 10;           // Adjusts vertical scale of spectrum plot


void qvisualize_transform_UART(_q * transformed)            // Sends plot to serial monitor
{
    int i;

    for(i = 0; i < UART_BUFFER_SIZE; i++)
    {
        memset((void *)UART_buffer[i], '_', ADC_BUFFER_SIZE/2+1);   // Set all pixels to just _
        UART_buffer[i][ADC_BUFFER_SIZE/2-1] = '\r';                 // add car. return and newline
        UART_buffer[i][ADC_BUFFER_SIZE/2] = '\n';
    }

    for (i = 0; i < ADC_BUFFER_SIZE/2-1; i++)
    {
        int j;
        // get next transformed value
        _q t_value = transformed[i+1];

        for (j = 0; j < UART_BUFFER_SIZE; j++)
        {
            // value to test if spectrum pixel should be populated
            int comp_value = ((spectrum_scale * (1 << GLOBAL_Q)  / UART_BUFFER_SIZE * j)-22000);

            if (t_value > comp_value)
            {
                // populate spectrum pixel
                UART_buffer[UART_BUFFER_SIZE-j-1][i] = '#';
            }
        }
    }

    UART_index = 0;
    UART_buffer_flag = 0;

    UCA1IE |= UCTXCPTIE;                // Turns on UART Tx IRQ
    UCA1IFG &= ~UCTXCPTIFG;             // Clear Tx complete flag
    UCA1TXBUF = UART_buffer[0][0];      // send first char of msg

    while (UART_buffer_flag == 0);      // wait until message is sent
}

int main(void)
{
    WDTCTL = WDTPW | WDTHOLD;           // stop watchdog timer

    //-- setup A1 for UART-------------------
    UCA1CTLW0 |= UCSWRST;               // put into SW reset

    UCA1CTLW0 |= UCSSEL__SMCLK;         // choose SMCLK (115200 baud)

    UCA1BRW = 8;                        // prescalar = 8
    UCA1MCTLW |= 0xD600;                // set modulation to to low frequency

    //-- setup ports------------------------
    P4SEL1 &= ~BIT3;                    // set P4.3 function to A1 UART Tx
    P4SEL0 |= BIT3;

    P1SEL1 |= BIT2;                     // configure P1.2 as Analog function
    P1SEL0 |= BIT2;

    P2REN |= BIT3;                      // Enable P2.3 resistors
    P4REN |= BIT1;                      // Enable P4.1 resistors

    P2OUT |= BIT3;                      // Set P2.3 resistors to pull up
    P4OUT |= BIT1;                      // Set P4.1 resistors to pull up


    P2IES |= BIT3;                      // Set P2.3 interrupt trigger to high to low
    P4IES |= BIT1;                      // Set P4.1 interrupt trigger to high to low

    PM5CTL0 &= ~LOCKLPM5;               // turn on gpio
    UCA1CTLW0 &= ~UCSWRST;              // take UART A1 out of SW reset

    UCB0CTLW0 &= ~UCSWRST;              // take B0 out of SW RST

    //-- Configure the ADC-----------------
    ADCCTL0 &= ~ADCSHT;                 // set conv. clock cycles = 16
    ADCCTL0 |= ADCSHT_2;
    ADCCTL0 |= ADCON;                   // turn on ADC

    ADCCTL1 |= ADCSSEL_2;               // choose SMCLK
    ADCCTL1 |= ADCSHP;                  // sample signal source = sampling timer

    ADCCTL2 &= ~ADCRES;                 // clear resolution
    ADCCTL2 |= ADCRES_2;                // 12-bit resolution

    ADCMCTL0 |= ADCINCH_2;              // ADC input = A2 (p1.2)

    //-- enable interrupts-----------------
    ADCIE |= ADCIE0;                    // Conversion complete IRQ

    P2IFG &= ~BIT3;                     // Clear P2.3 Interrupt flag
    P4IFG &= ~BIT1;                     // Clear P4.1 interrupt flag

    P2IE |= BIT3;                       // Enable P2.3 interrupts
    P4IE |= BIT1;                       // Enable P4.1 interrupts

    __enable_interrupt();               // enable maskables

    //-- local variables------------------
    int i;

    _q qre_buffer[ADC_BUFFER_SIZE];                                         // Stores real values of transform
    struct qcomplex qfreq_domain[ADC_BUFFER_SIZE];                          // Stores transformed result

    _q hamming_lookup[ADC_BUFFER_SIZE];                                     // Lookup table for hamming filter
    struct qcomplex euler_lookup[ADC_BUFFER_SIZE/2][ADC_BUFFER_SIZE+1];     // Lookup table for Euler's coef.

    populate_hamming_lookup(hamming_lookup, ADC_BUFFER_SIZE);               // Populate lookup tables
    populate_euler_lookup(euler_lookup, ADC_BUFFER_SIZE);


    while (1)
    {
        // clear fft buffers
        memset((void *)qre_buffer, 0, ADC_BUFFER_SIZE * sizeof(int));
        memset((void *)qfreq_domain, 0, ADC_BUFFER_SIZE * sizeof(struct qcomplex));

        for (i = 0; i < ADC_BUFFER_SIZE; i++)
        {
            // convert to volts and apply Hann window
            qre_buffer[i] = _Qmpy(hamming_lookup[i], _Q(ADC_buffer[i] * 3.4 / 4096 - 2.5));
        }

        ADC_index = 0;                          // reset index
        ADC_buffer_flag = 0;                    // reset buffer flag
        ADCCTL0 |= ADCENC | ADCSC;              // start conversion

        // Perform FFT
        qFFT_algorithm(euler_lookup, ADC_BUFFER_SIZE, qfreq_domain, qre_buffer, ADC_BUFFER_SIZE, 1);

        for (i = 0; i < ADC_BUFFER_SIZE; i++)
        {
            // Calculate amplitude of frequency on spectrum plot
            qre_buffer[i] = _Qlog(_Qmag(qfreq_domain[i].re, qfreq_domain[i].im));
        }

        qvisualize_transform_UART(qre_buffer);  // send the plot over serial

        while (ADC_buffer_flag == 0);           // wait for conversion to complete

    }

    return 0;
}

//--ISRs-----------------------------

#pragma vector = EUSCI_A1_VECTOR                        // UART Interrupt
__interrupt void ISR_EUSCI_A1(void)
{
    UCA1IFG &= ~UCTXCPTIFG;                             // clear Tx complete flag

    if (UART_index == (UART_BUFFER_SIZE)*(ADC_BUFFER_SIZE/2+1)-2)
    {
        UCA1IE &= ~UCTXCPTIE;                           // turn off Tx IRQ
        UART_buffer_flag = 1;
    }
    else
    {
        UART_index++;                                  // Increment position

        // Put next char into Tx buff
        UCA1TXBUF = UART_buffer[UART_index/(ADC_BUFFER_SIZE/2+1)][UART_index%(ADC_BUFFER_SIZE/2+1)];
    }
}

#pragma vector = ADC_VECTOR                             // ADC interrupt
__interrupt void ADC_ISR(void)
{
    ADC_buffer[ADC_index] = ADCMEM0;                    // Read ADC_Value
    ADC_index++;                                        // increment index

    if(ADC_index < ADC_BUFFER_SIZE)
    {
        ADCCTL0 |= ADCENC | ADCSC;                       // start next conversion
    }
    else
    {
        ADC_buffer_flag = 1;                            // mark that ADC sequence is complete
    }
}

#pragma vector = PORT2_VECTOR                           // Port 2 interrupt
__interrupt void ISR_Port2_S3(void)
{

    P2IFG &= ~BIT3;                                     // Clear P2.3 Interrupt flag


    spectrum_scale += spectrum_scale < 15 ? 1 : 0;      // Adjust scale

}

#pragma vector = PORT4_VECTOR                           // Port 1 interrupt
__interrupt void ISR_Port4_S1(void)
{

    P4IFG &= ~BIT1;                                     // Clear P4.1 Interrupt flag


    spectrum_scale -= spectrum_scale > 0 ? 1 : 0;       // Adjust scale
}

