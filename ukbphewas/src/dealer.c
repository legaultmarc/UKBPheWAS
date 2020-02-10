#include <string.h>
#include <stdio.h>
#include <zmq.h>

int main(void) {
    // The context
    void *context = zmq_ctx_new();

    // The required sockets (dealer and puller)
    void *dealer = zmq_socket(context, ZMQ_DEALER);
    void *puller = zmq_socket(context, ZMQ_PULL);

    // Binding the sockets
    zmq_bind(dealer, "tcp://*:0");
    zmq_bind(puller, "tcp://*:0");

    // Finding the endpoint for the dealer
    char dealer_endpoint[1024];
    size_t size_dealer_endpoint = sizeof(dealer_endpoint);
    zmq_getsockopt(dealer, ZMQ_LAST_ENDPOINT, dealer_endpoint, &size_dealer_endpoint);

    // Finding the endpoint for the puller
    char puller_endpoint[1024];
    size_t size_puller_endpoint = sizeof(puller_endpoint);
    zmq_getsockopt(puller, ZMQ_LAST_ENDPOINT, puller_endpoint, &size_puller_endpoint);

    printf("Starting to serve\nin=%s\nout_port=%s\n<>\n", puller_endpoint, dealer_endpoint);
    fflush(stdout);

    int nb_messages = 0;
    while (1) {
        zmq_msg_t message;
        while (1) {
            zmq_msg_init(&message);
            zmq_msg_recv(&message, puller, 0);

            int more = zmq_msg_more(&message);
            zmq_msg_send(&message, dealer, more? ZMQ_SNDMORE: 0);
            zmq_msg_close(&message);
            if (!more)
                break;
        }

        nb_messages++;
        if (nb_messages % 1000 == 0) printf("Sent %d messages\n", nb_messages);
    }

    zmq_close(dealer);
    zmq_close(puller);
    zmq_ctx_destroy(context);

    return 0;
}
